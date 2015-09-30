import sys
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time

# reference energy is from
# J.Sims, S.Hagstom, JCP *124*, 094101 (2006)
ref_E0 = -1.174475174
bond_length = 1.4
rmax = 60.0

# in reference paper
# order = 8
# rmax = 60
# num_knots =140
#lmax = 8

vars = [ (lmax, num, order)
         for lmax in [4]
         for num in [20, 60, 100, 140]
         for order in [6]]

acc = []
for (lmax, num, order) in vars:
    bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num))
    y_list = [CoupledY((L1, L2), L, 0)
              for L in range(lmax+1)
              for L1 in range(lmax+1)
              for L2 in range(lmax+1)
              if triangle_q(L1, L2, L) and (L1+L2)%2==0 and L1+L2<=lmax and L%2==0]
    t0 = time.clock()
    (h, s) = mat_h2(bond_length, bspline_set, y_list)
    t1 = time.clock()
    (es, vecs) = la.eigsh(h, 5, s, sigma=-3.0)
    t2 = time.clock()
    acc.append([("lmax", lmax),
                ("num", num),
                ("order", order),
                ("Nbais", h.shape[0]),
                ("Ndata", len(h.data)),
                ("t_mat", t1-t0),
                ("t_diag", t2-t1),
                ("E0", es[0] + 1.0/bond_length)])


for (name, val) in acc[0]:
    print name, " ",
print ""

for xs in acc:
    for (name, val) in xs:
        print val, " ",
    print ""
