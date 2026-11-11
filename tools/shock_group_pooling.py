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


def crossfit_pool_points(points):
    """Choose each peak using five other subgroups, then measure on the sixth.

    Density/width summaries remain unchanged. Profiles are independently
    centered to express training and evaluation peak coordinates consistently.
    This candidate needs separate coverage validation.
    """
    from shock_ensemble_estimator import centered_density,PROFILE_GRID
    from shock_crossfit_gradient import derivative,peak_position
    from shock_alignment import GRID
    if len(points)!=6:raise ValueError('Six independent subgroups required')
    result=pool_points(points)
    if set(result)!={'density','slopes','density_cubic','slope_cubic'}:
        raise ValueError('Both interpolation profiles required')
    profiles=np.asarray([p['density'] for p in points])
    alternatives=np.asarray([p['density_cubic'] for p in points])
    if profiles.shape!=(6,len(PROFILE_GRID)) or alternatives.shape!=profiles.shape:
        raise ValueError('Extended centered-profile support required')
    profiles=centered_density(profiles)
    alternatives=centered_density(alternatives,interpolation='cubic')
    slopes=[];alternate=[]
    for width in (1.,1.5,2.):
        values=[]
        for heldout in range(6):
            training=np.mean(np.delete(profiles,heldout,axis=0),axis=0)
            position=peak_position(GRID,training,width=width)
            x,d=derivative(profiles[heldout],width=width)
            slope=float(np.interp(position,x,d))
            if slope<=0:raise ValueError('Nonpositive heldout slope')
            values.append(slope)
            if width==1.5:
                training_cubic=np.mean(np.delete(alternatives,heldout,axis=0),axis=0)
                pos=peak_position(GRID,training_cubic,width=width)
                xx,dd=derivative(alternatives[heldout],width=width)
                s=float(np.interp(pos,xx,dd))
                if s<=0:raise ValueError('Nonpositive alternative slope')
                alternate.append(s)
        slopes.append(float(np.mean(values)))
    result['slopes']=slopes;result['slope_cubic']=float(np.mean(alternate))
    return result
