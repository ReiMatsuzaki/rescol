import sys
import commands
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time
import resource


# reference energy is from
# J.Sims, S.Hagstom, JCP *124*, 094101 (2006)
ref_E0 = -1.174475174
bond_length = 1.4
rmax = 60.0

vars = [ (lmax, num, order)
         for lmax in [2]
         for num in [20]
         for order in [3]]

t0 = time.clock()
acc = []
for (lmax, num, order) in vars:
    bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num))
    y_list = get_coupledY_set(lmax, 0, True, True)
    (h, s) = mat_h2(bond_length, bspline_set, y_list)
    (es, vecs) = la.eigsh(h, 1, s, sigma=-3.0)
    acc.append([("lmax", lmax),
                ("num", num),
                ("order", order),
                ("Nbasis", h.shape[0]),
                ("Ny2", len(y_list)),
                ("Ndata", len(h.data)),
                ("E0", es[0] + 1.0/bond_length)])

t1 = time.clock()

print "node: ", commands.getoutput("uname -n")
for (name, val) in acc[0]:
    print name, ": ", val
print "time: ", t1-t0
ru = resource.getrusage(resource.RUSAGE_SELF)
print "mem(Mb): ", ru.ru_maxrss/(1024.0)
