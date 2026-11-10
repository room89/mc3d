"""Paired controls for residual six-subgroup peak bias; never qualification.

Known fronts and peaks are available only for these synthetic controls.
All methods see identical random fixtures. The production estimator is unchanged.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_crossfit_gradient import derivative, peak_position
from shock_ensemble_estimator import group_point, centered_density
from shock_group_pooling import crossfit_pool_points
from diagnose_shock_ensemble_bias import quintic_derivative
from shock_validation import digest, save
from validate_shock_estimator import shape, truth


def slopes(profiles, position=None, quintic=False, width=1.5):
    profiles = np.asarray(profiles)
    if profiles.shape != (6,len(GRID)) or not np.isfinite(profiles).all(): raise ValueError('Six profiles required')
    values = []; peaks = []
    for heldout in range(6):
        training = np.delete(profiles,heldout,axis=0).mean(0)
        peak = peak_position(GRID,training,width=1.5) if position is None else position
        x, d = quintic_derivative(profiles[heldout],width=width) if quintic else derivative(profiles[heldout],width=width)
        values.append(float(np.interp(peak,x,d))); peaks.append(float(peak))
    if min(values)<=0: raise ValueError('Nonpositive diagnostic slope')
    return float(np.mean(values)), peaks


def trial(spec):
    kind, noise, repetition, windows, seed_start = spec
    position = math.log(math.log(2)) if kind=='gompertz' else 0.
    controls = {}; peak_errors = []; raw = []; cubic = []; oracles = []
    for group in range(12):
        points = []; oracle = []
        for subgroup in range(6):
            trials, exact = fixture(kind,seed_start+repetition*100+group*6+subgroup,0.,0.,noise=noise)
            trials = [r[:windows] for r in trials]; exact = exact[:,:windows]
            points.append(group_point(trials,center_profile=True,interpolation_sensitivity=True))
            oracle.append(exact.mean((0,1)))
        centered = centered_density([p['density'] for p in points])
        centered_cubic = centered_density([p['density_cubic'] for p in points],interpolation='cubic')
        uncentered = np.asarray([p['density'][20:-20] for p in points])
        production = crossfit_pool_points(points)
        measured, peaks = slopes(centered)
        if not np.isclose(measured,production['slopes'][1],rtol=1e-14,atol=0):
            raise ValueError('Control failed to reproduce production slope')
        values = dict(oracle_known_cubic=slopes(oracle,position)[0],
                      oracle_known_quintic=slopes(oracle,position,True)[0],
                      oracle_selected_cubic=slopes(oracle)[0],
                      refined_uncentered_cubic=slopes(uncentered)[0],
                      refined_centered_cubic=measured,
                      refined_centered_quintic_same_peaks=slopes(centered,quintic=True)[0],
                      refined_cubic_interpolation=production['slope_cubic'],
                      combined_quintic_15=slopes(centered_cubic,quintic=True)[0],
                      combined_quintic_20=slopes(centered_cubic,quintic=True,width=2.)[0],
                      combined_quintic_25=slopes(centered_cubic,quintic=True,width=2.5)[0])
        for name,value in values.items(): controls.setdefault(name,[]).append(value)
        peak_errors.extend(p-position for p in peaks)
        raw.extend(p['density'] for p in points); cubic.extend(p['density_cubic'] for p in points); oracles.extend(oracle)
    exact = truth(kind)[1]
    expected = shape(GRID-position,kind)
    density = dict(oracle=np.mean(oracles,axis=0),
                   refined_linear=centered_density(np.mean(raw,axis=0))[0],
                   refined_cubic=centered_density(np.mean(cubic,axis=0),interpolation='cubic')[0])
    return dict(repetition=repetition,relative_errors={k:float(1/np.mean(v)/exact-1) for k,v in controls.items()},
                peak_rmse=float(np.sqrt(np.mean(np.square(peak_errors)))),
                density_max_error={k:float(np.max(abs(v-expected))) for k,v in density.items()})


def cached(task):
    spec, folder, protocol = task; kind, noise, repetition, *_ = spec
    path = Path(folder)/f'{kind}-{noise}-{repetition}.json'
    if path.exists():
        saved = json.loads(path.read_text())
        if saved['protocol'] != protocol or saved['spec'] != list(spec): raise ValueError('Checkpoint mismatch')
        return saved['result']
    try: result = trial(spec)
    except ValueError as error: result = dict(repetition=repetition,error=str(error))
    save(path,dict(protocol=protocol,spec=spec,result=result)); return result


def summarize(rows):
    valid = [r for r in rows if 'error' not in r]; measures = {}
    def stat(values):
        a = np.asarray(values); return dict(mean=float(a.mean()),sd=float(a.std(ddof=1)),
                                           two_mc_se=float(2*a.std(ddof=1)/math.sqrt(len(a))))
    if len(valid)>=2:
        measures = {k:stat([r['relative_errors'][k] for r in valid]) for k in valid[0]['relative_errors']}
        contrasts = [('peak_selection','oracle_selected_cubic','oracle_known_cubic'),
                     ('alignment','refined_uncentered_cubic','oracle_selected_cubic'),
                     ('centering','refined_centered_cubic','refined_uncentered_cubic'),
                     ('derivative_order','refined_centered_quintic_same_peaks','refined_centered_cubic'),
                     ('interpolation','refined_cubic_interpolation','refined_centered_cubic')]
        for name in ('combined_quintic_15','combined_quintic_20','combined_quintic_25'):
            if name in valid[0]['relative_errors']: contrasts.append((name+'_minus_current',name,'refined_centered_cubic'))
        for name,left,right in contrasts:
            measures[name] = stat([r['relative_errors'][left]-r['relative_errors'][right] for r in valid])
    return dict(valid=len(valid),errors=[r for r in rows if 'error' in r],summary=measures,runs=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--replicates',type=int,default=8); p.add_argument('--jobs',type=int,default=8)
    p.add_argument('--windows',type=int,choices=(32,64),default=64); p.add_argument('--seed-start',type=int,default=100000000)
    p.add_argument('--noise-levels',type=float,nargs='+',default=[0.,.015])
    p.add_argument('--resume',action='store_true'); a=p.parse_args()
    if a.replicates<2 or a.jobs<1: p.error('At least two repetitions and one worker required')
    if not a.noise_levels or len(set(a.noise_levels))!=len(a.noise_levels) or any(not math.isfinite(n) or n<0 for n in a.noise_levels): p.error('Distinct finite nonnegative noise levels required')
    if a.output.exists(): raise FileExistsError(a.output)
    names=('diagnose_pooled_peak_bias.py','experiment_crossfit_gradient.py','shock_alignment.py',
           'shock_crossfit_gradient.py','shock_ensemble_estimator.py','shock_group_pooling.py',
           'diagnose_shock_ensemble_bias.py','shock_analysis.py','shock_mean_profile.py','validate_shock_estimator.py')
    protocol=dict(seed_start=a.seed_start,windows=a.windows,groups=12,subgroups=6,phi=0.,spatial=0.,noise_levels=a.noise_levels,
                  sources={n:digest(Path(__file__).with_name(n)) for n in names})
    parts=a.output.with_suffix('.parts'); parts.mkdir(parents=True,exist_ok=a.resume)
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for kind in ('logistic','gompertz'):
            for noise in a.noise_levels:
                specs=[(kind,noise,i,a.windows,a.seed_start) for i in range(a.replicates)]
                rows=list(pool.map(cached,[(s,str(parts),protocol) for s in specs]))
                cases.append(dict(shape=kind,noise=noise,**summarize(rows)))
                print(kind,noise,cases[-1]['valid'],cases[-1]['summary'],flush=True)
    if any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()): raise ValueError('Source changed')
    save(a.output,dict(status='diagnostic_only',protocol=protocol,cases=cases,
                      limitations=['Post-validation paired controls; not a fresh qualification.',
                                   'Known front/peak controls are unavailable for physical data.',
                                   'Independent group fixtures shared between paired controls; contrasts are conditional diagnostics.']))


if __name__=='__main__': main()
