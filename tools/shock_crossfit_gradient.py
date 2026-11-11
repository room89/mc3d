"""Experimental independent peak selection; never maximize on the held-out seed."""
import numpy as np
from shock_alignment import GRID


def derivative(profile, x=GRID, width=1.):
    x=np.asarray(x,float);y=np.asarray(profile,float)
    if len(x)<5 or len(x)!=len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite profile required')
    dx=x[1]-x[0]
    if dx<=0 or not np.allclose(np.diff(x),dx):raise ValueError('Uniform increasing grid required')
    radius=int(round(width/(2*dx)))
    if radius<2 or 2*radius+1>len(x):raise ValueError('Unresolved derivative')
    t=np.arange(-radius,radius+1)*dx/width
    weights=np.linalg.pinv(np.column_stack([t**k for k in range(4)]))[1]/width
    return x[radius:-radius],np.convolve(y,weights[::-1],mode='valid')


def peak_position(template_x,template_density,width=1.):
    x,d=derivative(template_density,template_x,width)
    index=int(d.argmax())
    if index in (0,len(d)-1) or d[index]<=0:raise ValueError('Invalid template peak')
    return float(x[index])


def estimate_gradient(data,diagnostics,width=1.):
    data=np.asarray(data,float)
    if data.ndim!=3 or data.shape[0]!=3 or data.shape[2]!=len(GRID) or not np.isfinite(data).all() or len(diagnostics)!=3:
        raise ValueError('Three aligned seeds and their templates required')
    slopes=[];details=[]
    for seed,(run,info) in enumerate(zip(data,diagnostics)):
        if info['training_seeds']!=[s for s in range(3) if s!=seed]:raise ValueError('Template must exclude held-out seed')
        position=peak_position(info['template_x'],info['template_density'],width)
        x,d=derivative(run.mean(axis=0),GRID,width)
        if not x[0]<=position<=x[-1]:raise ValueError('Peak outside held-out coverage')
        slope=float(np.interp(position,x,d))
        if slope<=0:raise ValueError('Nonpositive held-out derivative')
        slopes.append(slope);details.append(dict(heldout_seed=seed,peak_position=position,slope=slope))
    return dict(gradient_thickness=float(1/np.mean(slopes)),folds=details)
