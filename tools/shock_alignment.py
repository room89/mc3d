"""Cross-fitted translation-only alignment. All lengths are in upstream lambda."""
from dataclasses import dataclass
import numpy as np
from shock_analysis import monotone,crossing

GRID=np.linspace(-10,10,401)

@dataclass(frozen=True)
class PreparedWindow:
    x: np.ndarray
    y: np.ndarray
    crossing50: float


def prepare_trials(trials):
    cache={}
    def prepare(w):
        if isinstance(w,PreparedWindow):return w
        key=id(w)
        if key not in cache:
            x,y,c=initial(w)
            x=x.copy();y=y.copy();x.setflags(write=False);y.setflags(write=False)
            cache[key]=PreparedWindow(x,y,c)
        return cache[key]
    return [[prepare(w) for w in run] for run in trials]


def initial(window):
    if isinstance(window,PreparedWindow):return window.x,window.y,window.crossing50
    y=np.asarray(window['density_normalized'],float)
    x=np.asarray(window['x'] if 'x' in window else window['x_aligned'],float)
    if 'x' not in window:x=x+window.get('front',0.)
    if len(x)!=len(y) or len(x)<4 or not np.isfinite(x).all() or not np.isfinite(y).all() or not np.all(np.diff(x)>0):raise ValueError('Invalid profile')
    c=crossing(x,monotone(y),.5)
    return x,y,c

def interp(x,y,points):
    if points[0]<x[0] or points[-1]>x[-1]:raise ValueError('Insufficient profile coverage')
    return np.interp(points,x,y)

def smooth_template(rows):
    offsets=np.arange(-15,16)*.05;kernel=np.exp(-.5*(offsets/.25)**2);kernel/=kernel.sum()
    mean=np.mean(rows,axis=0)
    smooth=np.convolve(np.pad(mean,15,mode='edge'),kernel,mode='valid')
    smooth=np.asarray(monotone(smooth));center=crossing(GRID,smooth,.5)
    return GRID-center,smooth


def template(trials,heldout,refine_passes=0):
    if refine_passes not in (0,1):raise ValueError('Supported refinement passes: zero or one')
    rows=[];training=[]
    for seed,run in enumerate(trials):
        if seed==heldout:continue
        for w in run:
            x,y,c=initial(w);rows.append(interp(x,y,GRID+c));training.append(w)
    tx,ty=smooth_template(rows)
    if refine_passes:
        rows=[]
        for w,(front,_) in zip(training,shift_many(training,tx,ty)):
            x,y,_=initial(w)
            rows.append(interp(x,y,GRID+front))
        tx,ty=smooth_template(rows)
    return tx,ty

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

def shift_many(windows,tx,ty):
    """Vectorized equivalent of shift on a shared source grid; scalar fallback."""
    raw=[initial(w) for w in windows]
    if not raw:return []
    x=raw[0][0]
    if any(len(v[0])!=len(x) or not np.allclose(v[0],x,rtol=1e-14,atol=1e-12) for v in raw):
        return [shift(w,tx,ty) for w in windows]
    y=np.array([v[1] for v in raw]);c=np.array([v[2] for v in raw])
    mask=(ty>=.1)&(ty<=.9);u=tx[mask];target=ty[mask]
    if len(u)<4:raise ValueError('Unresolved template')
    def loss(delta):
        # delta shape (windows, candidate shifts).
        points=u[None,None,:]+c[:,None,None]+delta[:,:,None]
        if points.min()<x[0] or points.max()>x[-1]:raise ValueError('Insufficient profile coverage')
        j=np.clip(np.searchsorted(x,points,side='right')-1,0,len(x)-2)
        fraction=(points-x[j])/(x[j+1]-x[j])
        row=np.arange(len(y))[:,None,None]
        values=y[row,j]+fraction*(y[row,j+1]-y[row,j])
        residual=abs(values-target)
        return np.mean(np.where(residual<=.02,.5*residual**2,.02*(residual-.01)),axis=-1)
    candidates=np.linspace(-1,1,41)
    cost=loss(np.broadcast_to(candidates,(len(y),41)));indices=cost.argmin(axis=1)
    if np.any((indices==0)|(indices==40)):raise ValueError('Boundary minimum')
    for row,i in zip(cost,indices):
        for j in range(1,40):
            if abs(candidates[j]-candidates[i])>=.25 and row[j]<=row[j-1] and row[j]<=row[j+1] and row[j]-row[i]<=max(1e-8,.01*row[i]):raise ValueError('Ambiguous minimum')
    lo=candidates[indices-1];hi=candidates[indices+1];ratio=(5**.5-1)/2
    while np.max(hi-lo)>1e-4:
        a=hi-ratio*(hi-lo);b=lo+ratio*(hi-lo)
        left=loss(a[:,None])[:,0]<loss(b[:,None])[:,0]
        hi=np.where(left,b,hi);lo=np.where(left,lo,a)
    delta=(lo+hi)/2;cost=loss(delta[:,None])[:,0]
    return [(float(v+d),dict(initial_front=float(v),correction=float(d),front=float(v+d),loss=float(e),status='ok')) for v,d,e in zip(c,delta,cost)]


def align(trials,grid=GRID,refine_passes=0):
    if len(trials)!=3 or len({len(r) for r in trials})!=1:raise ValueError('Three equal-length seeds required')
    trials=prepare_trials(trials)
    result=[];diagnostics=[]
    for seed,run in enumerate(trials):
        tx,ty=template(trials,seed,refine_passes);rows=[];windows=[]
        for w,(front,detail) in zip(run,shift_many(run,tx,ty)):
            x,y,_=initial(w)
            rows.append(interp(x,y,grid+front));windows.append(detail)
        result.append(rows);diagnostics.append(dict(training_seeds=[i for i in range(3) if i!=seed],refine_passes=refine_passes,template_x=tx.tolist(),template_density=ty.tolist(),windows=windows))
    return np.asarray(result),diagnostics
