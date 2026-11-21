"""Independent SPARTA VHS equilibrium event-rate gate."""
import argparse,json,math,subprocess,hashlib
from pathlib import Path
from argon_vhs import fit,MASS,KB
from normal_shock import downstream
from shock_validation import save

def rate(t,gas):
    w=gas['vhs_omega'];d=gas['vhs_diameter_m'];p=2-2*w
    prefactor=math.pi*d*d*(4*KB*300/MASS)**(w-.5)/math.gamma(2.5-w)
    return prefactor*(4*KB*t/MASS)**(p/2)*math.gamma((p+3)/2)/math.gamma(1.5)

def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    gas,_=fit();results=[]
    for t in (300.,downstream(3.8,1e20,300)['temperature_k']):
        for seed in (123,456,789):
            folder=a.output_dir/f'{t:.3f}-{seed}';folder.mkdir();N=5000;n=1e20;dt=.02/(n*rate(t,gas))
            (folder/'argon.species').write_text(f'Ar 39.948 {MASS:.17g} 0 0 0 0 0 1 0\n')
            (folder/'argon.vss').write_text(f"Ar {gas['vhs_diameter_m']:.17g} {gas['vhs_omega']:.17g} 300 1\n")
            text=f'''seed {seed}
dimension 3
boundary p p p
create_box 0 .01 0 .01 0 .01
create_grid 1 1 1
species argon.species Ar
mixture gas Ar nrho {n} temp {t} vstream 0 0 0
global fnum {n*1e-6/N}
create_particles gas n {N}
collide vss all argon.vss
compute temp temp
stats 1
stats_style step np ncoll c_temp
timestep {dt}
run 1000
'''
            (folder/'in.rate').write_text(text)
            with (folder/'run.log').open('w') as log:subprocess.run([str(a.binary.resolve()),'-in','in.rate'],cwd=folder,stdout=log,stderr=log,check=True,timeout=120)
            rows=[]
            for line in (folder/'run.log').read_text().splitlines():
                fields=line.split()
                if len(fields)==4 and fields[0].isdigit():
                    try:rows.append([float(x) for x in fields])
                    except ValueError:pass
            rows=[r for r in rows if r[0]>0]
            if len(rows)!=1000:raise ValueError('Incomplete SPARTA statistics')
            observed=sum(r[2] for r in rows);actual_t=sum(r[3] for r in rows)/len(rows)
            expected=N*(N-1)/2*(n/N)*rate(actual_t,gas)*dt*len(rows)
            # Six-sigma counting plus finite-ensemble fluctuation allowance.
            tolerance=6*math.sqrt(expected)+expected*6*math.sqrt(2/(3*N))
            results.append(dict(temperature_k=t,seed=seed,observed=observed,expected=expected,tolerance=tolerance,passed=abs(observed-expected)<=tolerance))
    save(a.output_dir/'report.json',dict(binary_sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),runs=results))
    return 0 if all(r['passed'] for r in results) else 1
if __name__=='__main__':raise SystemExit(main())
