"""Diagnostic comparison of secondary experimental reproductions; no acceptance gate."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from shock_validation import digest,save
from shock_analysis import monotone


def crossings(points,coordinate,scale=1.):
    if not math.isfinite(scale) or scale<=0:raise ValueError('Positive finite scale required')
    x=np.asarray([p[coordinate] for p in points],float)*scale
    y=np.asarray([p['rho_normalized'] for p in points],float)
    if len(x)<3 or not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(np.diff(x)<=0):
        raise ValueError('Finite increasing coordinates required')
    y=monotone(y)
    if y[0]>.1 or y[-1]<.9:raise ValueError('Missing 10–90% support')
    return np.interp([.1,.5,.9],y,x)


def compare_sources(normalized,dimensional,reference_mean_free_path):
    a=crossings(normalized,'x_paper',reference_mean_free_path)
    b=crossings(dimensional,'x_m')
    wa=float(a[2]-a[0]);wb=float(b[2]-b[0])
    return dict(status='secondary_source_diagnostic_only',reference_mean_free_path_m=reference_mean_free_path,
        normalized_source_width_m=wa,dimensional_source_width_m=wb,
        dimensional_to_normalized_width_relative_difference=wb/wa-1,
        midpoint_translation_m=float(b[1]-a[1]),
        implied_length_from_width_ratio_m=reference_mean_free_path*wb/wa,
        centered_crossing_differences_m=((b-b[1])-(a-a[1])).tolist(),
        method='Monotone regression of markers then linear 10%, 50%, 90% crossings; fixed supplied length scale. Translation removed only for centered crossing diagnostic.',
        limitations=['Neither reproduction is primary experimental data.',
                    'Implied length from width ratio is diagnostic only, never an accepted fitted spatial scale.',
                    'Source discrepancy is not an experimental confidence interval or a justified acceptance tolerance.',
                    'Original coordinate convention and digitization uncertainty require verification.'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--normalized',type=Path,required=True);p.add_argument('--dimensional',type=Path,required=True)
    p.add_argument('--reference-mean-free-path-m',type=float,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    result=compare_sources(json.loads(a.normalized.read_text())['points'],json.loads(a.dimensional.read_text())['points'],a.reference_mean_free_path_m)
    result['source_hashes']={'normalized':digest(a.normalized),'dimensional':digest(a.dimensional)}
    result['tool_sha256']=digest(Path(__file__))
    save(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
