from sympy import *

"""
To derive exact value for dipole integral in test_oce1.cpp
"""

h1_2 = Rational(1,2)

r = Symbol('r')
r1s = 2 * r * exp(-r)
r2p = 1/(2*sqrt(6)) * r**2 * exp(-r/2)
print r2p
print "<R1s|R1s> = ", integrate(r1s*r1s, (r, 0, oo))
print "<R2p|R2p> = ", integrate(r2p*r2p, (r, 0, oo))
print "<R2p|r|R1s> = ", integrate(r2p*r*r1s, (r, 0, oo))

""" output:
<R1s|R1s> =  1
<R2p|R2p> =  1
<R2p|r|R1s> =  128*sqrt(6)/243
"""


