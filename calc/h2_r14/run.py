import sys
sys.path.append("../../src")
from rescol import *
import scipy.sparse.linalg as la
import time

# reference energy is from
# J.Sims, S.Hagstom, JCP *124*, 094101 (2006)
ref_E0 = -1.174475174
bond_length = 1.4


# in reference paper
# order = 8
# rmax = 60
# num_knots =140
#lmax = 8

order = 4
rmax = 60.0
num_knots = 21
lmax = 2

bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num_knots))
y_list = [CoupledY((L1, L2), L, 0)
          for L in range(2*lmax+1)
          for L1 in range(lmax+1)
          for L2 in range(lmax+1)
          if triangle_q(L1, L2, L) and (L1+L2)%2==0]

((h, s), calc_log) = mat_h2(bond_length, bspline_set, y_list)
t0 = time.clock()
(es, vecs) = la.eigsh(h, 5, s, sigma=-3.0)
t1 = time.clock()

calc_log = [("E0", es[0] + 1.0/bond_length)] + calc_log
calc_log.append(("diagonalize", t1-t0))
for (key, val) in calc_log:
    print key, ": ", val


