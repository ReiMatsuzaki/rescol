from sympy import *

"""
To derive exact value for dipole integral in test_oce1.cpp
"""

h1_2 = Rational(1,2)

r = Symbol('r')
r1s = 2 * r * exp(-r)
mu_v_r1s = 2 * r * exp(-r)
r2p = 1/(2*sqrt(6)) * r**2 * exp(-r/2)
mu_v_r2pplus = -1/(4*sqrt(6)) * r**2 * exp(-r/2)
r3d = 4/(81*sqrt(30)) * r**3 * exp(-r/3)
D = lambda f: diff(f, r)
print r2p
print "<R1s|R1s> = ", integrate(r1s*r1s, (r, 0, oo))
print "<R2p|R2p> = ", integrate(r2p*r2p, (r, 0, oo))
print "<R2p|r|R1s> = ", integrate(r2p*r*r1s, (r, 0, oo))
print "<R2p|dr|R1s> = " , integrate(r2p*D(r1s), (r,0,oo))
print "<R2p|mu_v_R1s> = " , integrate(r2p*mu_v_r1s, (r,0,oo))
print "<R3d|mu_v_R2p> = " , integrate(r3d*mu_v_r2pplus, (r,0,oo))

""" output:
<R1s|R1s> =  1
<R2p|R2p> =  1
<R2p|r|R1s> =  128*sqrt(6)/243
<R2p|dr|R1s> =  -8*sqrt(6)/81
<R2p|mu_v_R1s> =  16*sqrt(6)/81
<R3d|mu_v_R2p> =  -2304*sqrt(5)/15625
"""


