"""Offline argon fit: Sotiriadou et al. (2025), Eq. 2, Table 3.
DOI: 10.1007/s10765-025-03603-8; dilute limit valid 84–10000 K.
Run: python3 tools/argon_vhs.py --output DIR
"""
import argparse
import json
import math
from pathlib import Path

COEFFICIENTS = (.8395115, -.1062564, .01065796, .01879809,
                -.008881774, -.00009613779, .001404406, -.0004321739,
                -.00002544782, .00004398471, -.000009997908, .0000007753453)
MASS = 39.948 * 1.66053906892e-27
KB = 1.380649e-23


def viscosity(t):
    if not 84 <= t <= 10000:
        raise ValueError('Correlation valid only at 84–10000 K')
    x = math.log(t / 298.15)
    return 22.5666e-6 * math.exp(sum(a*x**i for i, a in enumerate(COEFFICIENTS, 1)))


def fit():
    xs = [math.log(2000/300)*i/31 for i in range(32)]
    ys = [math.log(viscosity(300*math.exp(x))) for x in xs]
    xm, ym = sum(xs)/32, sum(ys)/32
    omega = sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/sum((x-xm)**2 for x in xs)
    mu = math.exp(ym-omega*xm)
    diameter = math.sqrt(15*math.sqrt(MASS*KB*300/math.pi)/(2*(5-2*omega)*(7-2*omega)*mu))
    rows = []
    for i in range(1001):
        t = 300*math.exp(math.log(2000/300)*i/1000)
        ref = viscosity(t)
        value = mu*(t/300)**omega
        rows.append(dict(temperature_k=t, reference_pa_s=ref, vhs_pa_s=value, relative_error=value/ref-1))
    worst = max(rows, key=lambda r: abs(r['relative_error']))
    config = dict(collision_model='vhs', molecular_mass_kg=MASS, vhs_diameter_m=diameter,
                  vhs_reference_temperature_k=300, vhs_omega=omega, temperature_unit_k=300)
    report = dict(source='https://doi.org/10.1007/s10765-025-03603-8', coefficients=COEFFICIENTS,
                  mu_reference_pa_s=mu, max_relative_error=abs(worst['relative_error']),
                  max_error_temperature_k=worst['temperature_k'],
                  rms_relative_error=math.sqrt(sum(r['relative_error']**2 for r in rows)/len(rows)),
                  accepted=abs(worst['relative_error']) <= .05 and .5 <= omega <= 1,
                  fit_samples=[dict(temperature_k=300*math.exp(x), viscosity_pa_s=math.exp(y)) for x,y in zip(xs,ys)],
                  validation_samples=rows)
    return config, report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config, report = fit()
    args.output.mkdir(parents=True, exist_ok=True)
    for name, data in [('argon_vhs.json', config), ('argon_vhs_report.json', report)]:
        (args.output/name).write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({**config, **{k:v for k,v in report.items() if 'samples' not in k and k != 'coefficients'}}, indent=2))
    raise SystemExit(0 if report['accepted'] else 1)
