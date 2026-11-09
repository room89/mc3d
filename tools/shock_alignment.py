"""Cross-fitted translation-only alignment. All lengths are in upstream lambda."""
import numpy as np
from shock_analysis import monotone,crossing

GRID=np.linspace(-10,10,401)

def initial(window):
    y=np.asarray(window['density_normalized'],float)
    x=np.asarray(window.get('x',window['x_aligned']),float)
    if 'x' not in window:x=x+window.get('front',0.)
    if len(x)!=len(y) or len(x)<4 or not np.isfinite(x).all() or not np.isfinite(y).all() or not np.all(np.diff(x)>0):raise ValueError('Invalid profile')
    c=crossing(x,monotone(y),.5)
    return x,y,c

def interp(x,y,points):
    if points[0]<x[0] or points[-1]>x[-1]:raise ValueError('Insufficient profile coverage')
    return np.interp(points,x,y)

def template(trials,heldout):
    rows=[]
    for seed,run in enumerate(trials):
        if seed==heldout:continue
        for w in run:
            x,y,c=initial(w);rows.append(interp(x,y,GRID+c))
    offsets=np.arange(-15,16)*.05;kernel=np.exp(-.5*(offsets/.25)**2);kernel/=kernel.sum()
    mean=np.mean(rows,axis=0)
    smooth=np.convolve(np.pad(mean,15,mode='edge'),kernel,mode='valid')
    smooth=np.asarray(monotone(smooth));center=crossing(GRID,smooth,.5)
    return GRID-center,smooth

def shift(window,tx,ty):
    x,y,c=initial(window);mask=(ty>=.1)&(ty<=.9);u=tx[mask];target=ty[mask]
    if len(u)<4:raise ValueError('Unresolved template')
    def loss(d):
        residual=abs(interp(x,y,u+c+d)-target)
        return float(np.mean(np.where(residual<=.02,.5*residual**2,.02*(residual-.01))))
    candidates=np.linspace(-1,1,41);cost=np.array([loss(d) for d in candidates]);i=int(cost.argmin())
    if i in (0,40):raise ValueError('Boundary minimum')
    for j in range(1,40):
        if abs(candidates[j]-candidates[i])>=.25 and cost[j]<=cost[j-1] and cost[j]<=cost[j+1] and cost[j]-cost[i]<=max(1e-8,.01*cost[i]):raise ValueError('Ambiguous minimum')
    lo,hi=candidates[i-1],candidates[i+1];ratio=(5**.5-1)/2
    while hi-lo>1e-4:
        a=hi-ratio*(hi-lo);b=lo+ratio*(hi-lo)
        if loss(a)<loss(b):hi=b
        else:lo=a
    delta=(lo+hi)/2
    return c+delta,dict(initial_front=c,correction=delta,front=c+delta,loss=loss(delta),status='ok')

def align(trials,grid=GRID):
    if len(trials)!=3 or len({len(r) for r in trials})!=1:raise ValueError('Three equal-length seeds required')
    result=[];diagnostics=[]
    for seed,run in enumerate(trials):
        tx,ty=template(trials,seed);rows=[];windows=[]
        for w in run:
            front,detail=shift(w,tx,ty);x,y,_=initial(w)
            rows.append(interp(x,y,grid+front));windows.append(detail)
        result.append(rows);diagnostics.append(dict(training_seeds=[i for i in range(3) if i!=seed],template_x=tx.tolist(),template_density=ty.tolist(),windows=windows))
    return np.asarray(result),diagnostics
