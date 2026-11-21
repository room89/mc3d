import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from analytical_reference import cube_flux
class AnalyticalReferenceTest(unittest.TestCase):
    def test_reference(self):
        f,q=cube_flux()
        self.assertAlmostEqual(f[0],2.845842437,places=8)
        self.assertAlmostEqual(q,8.808258290,places=8)
        self.assertEqual(f[1:],[0.,0.])
    def test_equilibrium(self):
        f,q=cube_flux(speed=0)
        for x in f:self.assertAlmostEqual(x,0,places=14)
        self.assertAlmostEqual(q,0,places=14)
    def test_reversed_flow(self):
        f,q=cube_flux();r,h=cube_flux(speed=-5*math.sqrt(2))
        self.assertAlmostEqual(f[0],-r[0]);self.assertAlmostEqual(q,h)
    def test_independent_quadrature(self):
        # Integrate the incoming normal-speed density directly (Simpson rule).
        def integral(u,power):
            n=20000; upper=abs(u)+12; h=upper/n
            def f(v):return v**power*math.exp(-(v-u)**2/2)/math.sqrt(2*math.pi)
            return h/3*(f(0)+f(upper)+sum((4 if i%2 else 2)*f(i*h) for i in range(1,n)))
        U=5*math.sqrt(2); fx=0; heat=0
        for axis in range(3):
            for sign in (-1,1):
                u=-U*sign if axis==0 else 0
                a,b,c=(integral(u,k) for k in (1,2,3))
                fx+=.04*(-sign*(b+a*math.sqrt(math.pi/2)) if axis==0 else a*U)
                heat+=.04*(.5*(c+a*((0 if axis==0 else U*U)+2))-2*a)
        f,q=cube_flux();self.assertAlmostEqual(f[0],fx,places=9);self.assertAlmostEqual(q,heat,places=9)
