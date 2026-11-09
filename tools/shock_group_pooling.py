"""Pool equally sized independent three-seed subgroup sufficient statistics.

Average densities and slopes before nonlinear width/inverse-slope inference.
This is a candidate protocol; it requires its own coverage validation.
"""
import numpy as np


def pool_points(points):
    if not points:raise ValueError('At least one subgroup required')
    keys=set(points[0])
    if keys not in ({'density','slopes'},{'density','slopes','density_cubic','slope_cubic'}):
        raise ValueError('Unexpected subgroup statistics')
    if any(set(p)!=keys for p in points):raise ValueError('Mixed subgroup protocols')
    result={}
    for key in sorted(keys):
        arrays=[np.asarray(p[key],dtype=float) for p in points]
        if any(a.shape!=arrays[0].shape or not np.isfinite(a).all() for a in arrays):
            raise ValueError('Inconsistent or nonfinite subgroup statistics')
        result[key]=np.mean(arrays,axis=0).tolist()
    return result
