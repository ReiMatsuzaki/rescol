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

vars = [ (lmax, num, order)
         for lmax in [6]
         for num in [100]
         for order in [6]]

acc = []
for (lmax, num, order) in vars:
    bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num))
    y_list = get_coupledY_set(lmax, 0, True, True)
    t0 = time.clock()
    (h, s) = mat_h2(bond_length, bspline_set, y_list)
    t1 = time.clock()
    (es, vecs) = la.eigsh(h, 5, s, sigma=-3.0)
    t2 = time.clock()
    acc.append([("lmax", lmax),
                ("num", num),
                ("order", order),
                ("Nbasis", h.shape[0]),
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
