"""Independent Maxwell half-space flux reference, rho=1, particle mass=1."""
import math

def cube_flux(speed=5*math.sqrt(2), temperature=1., wall_temperature=1., side=.2):
    force=[0.,0.,0.]; heat=0.
    for axis in range(3):
        for sign in (-1,1):
            u=-speed*sign if axis==0 else 0.
            tangential=[speed,0.,0.]; tangential[axis]=0.
            phi=math.exp(-u*u/(2*temperature))/math.sqrt(2*math.pi)
            cdf=.5*math.erfc(-u/math.sqrt(2*temperature))
            j0=math.sqrt(temperature)*phi+u*cdf
            j1=(u*u+temperature)*cdf+u*math.sqrt(temperature)*phi
            j2=(u**3+3*u*temperature)*cdf+(u*u+2*temperature)*math.sqrt(temperature)*phi
            for k in range(3):
                force[k]+=side**2*(j0*tangential[k]-(sign if k==axis else 0)*(j1+j0*math.sqrt(math.pi*wall_temperature/2)))
            heat+=side**2*(.5*(j2+j0*(sum(v*v for v in tangential)+2*temperature))-2*wall_temperature*j0)
    return force,heat
