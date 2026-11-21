"""Pure-Python profile analysis; no physical acceptance from execution alone."""
import math,statistics

def monotone(values):
    blocks=[]
    for value in values:
        if not math.isfinite(value):raise ValueError('Nonfinite density')
        blocks.append([value,1])
        while len(blocks)>1 and blocks[-2][0]/blocks[-2][1]>blocks[-1][0]/blocks[-1][1]:
            b=blocks.pop();blocks[-1][0]+=b[0];blocks[-1][1]+=b[1]
    return [total/n for total,n in blocks for _ in range(n)]

def crossing(x,y,target):
    if len(x)!=len(y) or len(x)<2 or any(b<=a for a,b in zip(x,x[1:])):raise ValueError('Invalid coordinates')
    for i in range(len(x)-1):
        if y[i]<=target<y[i+1]:return x[i]+(target-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])
    raise ValueError('Profile does not contain requested crossing')

def solve(a,b):
    a=[list(row)+[v] for row,v in zip(a,b)];n=len(b)
    for k in range(n):
        pivot=max(range(k,n),key=lambda i:abs(a[i][k]));a[k],a[pivot]=a[pivot],a[k]
        scale=a[k][k]
        if abs(scale)<1e-15:raise ValueError('Singular local fit')
        a[k]=[v/scale for v in a[k]]
        for i in range(n):
            if i!=k:
                scale=a[i][k];a[i]=[v-scale*w for v,w in zip(a[i],a[k])]
    return [row[-1] for row in a]

def gradient(x,y,width):
    if width<=0:raise ValueError('Invalid fit width')
    result=[]
    for center in x:
        pairs=[((a-center)/width,b) for a,b in zip(x,y) if abs(a-center)<=width/2]
        if len(pairs)<4:continue
        a=[[sum(t**(i+j) for t,_ in pairs) for j in range(4)] for i in range(4)]
        b=[sum(v*t**i for t,v in pairs) for i in range(4)]
        result.append(solve(a,b)[1]/width)
    if not result or max(result)<=0:raise ValueError('Unresolved gradient')
    return max(result)

def metrics(x,density,n1,n2,mean_free_path):
    if n2<=n1 or mean_free_path<=0:raise ValueError('Invalid shock reference')
    y=[(v-n1)/(n2-n1) for v in density];m=monotone(y)
    front=crossing(x,m,.5);width=crossing(x,m,.9)-crossing(x,m,.1)
    thick=1/gradient(x,y,mean_free_path)
    alternative=1/gradient(x,y,mean_free_path/2)
    return dict(front=front,width_10_90=width,gradient_thickness=thick,processing_uncertainty=abs(thick-alternative),x_aligned=[v-front for v in x],density_normalized=y)

def interval(values):
    if len(values)<3 or not all(math.isfinite(v) for v in values):raise ValueError('Need at least three finite independent samples')
    # Conservative t critical value: df=2 for all n>=3.
    return statistics.mean(values),4.303*statistics.stdev(values)/math.sqrt(len(values))

def equivalence(difference,uncertainty,tolerance):
    if not all(math.isfinite(v) for v in (difference,uncertainty,tolerance)) or uncertainty<0 or tolerance<=0:raise ValueError('Invalid comparison')
    return uncertainty<=tolerance/2 and abs(difference)+uncertainty<=tolerance

def pooled(rows,weight,duration,volume):
    if weight<=0 or duration<=0 or volume<=0:raise ValueError('Invalid normalization')
    a0=sum(float(r['A0']) for r in rows)
    if a0<=0:raise ValueError('Empty profile bin')
    axes='xyz';u=[sum(float(r['A1'+a])+float(r['origin_'+a])*float(r['A0']) for r in rows)/a0 for a in axes]
    tensor=[[0.]*3 for _ in range(3)];q=[0.]*3
    for r in rows:
        n=float(r['A0']);a1=[float(r['A1'+a]) for a in axes];a3=[float(r['A3'+a]) for a in axes]
        a2=[[float(r['A2'+''.join(sorted(a+b))]) for b in axes] for a in axes]
        d=[float(r['origin_'+a])-u[i] for i,a in enumerate(axes)];d2=sum(v*v for v in d);dot=sum(d[i]*a1[i] for i in range(3));tr=sum(a2[i][i] for i in range(3))
        for i in range(3):
            q[i]+=a3[i]+d[i]*tr+2*sum(a2[i][j]*d[j] for j in range(3))+2*d[i]*dot+d2*a1[i]+n*d2*d[i]
            for j in range(3):tensor[i][j]+=a2[i][j]+d[i]*a1[j]+d[j]*a1[i]+n*d[i]*d[j]
    scale=weight/(duration*volume)
    return dict(density=scale*a0,velocity=u,temperature=sum(tensor[i][i] for i in range(3))/(3*a0),pressure=[[v*scale for v in row] for row in tensor],heat=[v*scale/2 for v in q])
