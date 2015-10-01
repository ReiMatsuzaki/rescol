import sys
import commands
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time
import resource
import traceback


# reference energy is from
# J.Sims, S.Hagstom, JCP *124*, 094101 (2006)
ref_E0 = -1.174475174
bond_length = 1.4
rmax = 60.0

vars = [ (lmax, num, order)
         for lmax in [4]
         for num in [100]
         for order in [6]]

t0 = time.clock()
print "node: ", commands.getoutput("uname -n")
for (lmax, num, order) in vars:
    print "lmax: ", lmax
    print "num: ", num
    print "order: ", order

    bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num))
    y_list = get_coupledY_set(lmax, 0, True, True)
    try:
        (h, s) = mat_h2(bond_length, bspline_set, y_list)
    except:
        print "Failed at building matrix"
        print "--------------------------------"
        print traceback.format_exc(sys.exc_info()[2])
        print "--------------------------------"

    print "NBasis: ", h.shape[0]
    print "Ny2: ", len(y_list)
    print "NData: {0:e}".format(len(h.data))
    write_coo_mat(h, "h.out")
    write_coo_mat(s, "s.out")

    """
    try:
        (es, vecs) = la.eigsh(h, 1, s, sigma=-3.0)
    except:
        print "Failed at diagonalizing matrix"
        print "--------------------------------"
        print traceback.format_exc(sys.exc_info()[2])
        print "--------------------------------"
    print "E0: ", es[0]+1.0/bond_length
    """

t1 = time.clock()
print "time: ", t1-t0
ru = resource.getrusage(resource.RUSAGE_SELF)
print "mem(Mb): ", ru.ru_maxrss/(1024.0)
