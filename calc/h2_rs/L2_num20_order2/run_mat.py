import sys
import os
import commands
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time
import resource
import traceback

bond_length_list = [i*0.1+0.3 for i in range(40)]
rmax = 20.0

vars = [ (lmax, num, order)
         for lmax in [2]
         for num in [20]
         for order in [2]]

os.chdir("calc")
for bond_length in bond_length_list:
    dir_name = "{0}".format(round(bond_length*10))
    print dir_name
    if not os.path.exists(dir_name):
        os.mkdir(dir_name) 
    t0 = time.clock()
    os.chdir(dir_name)
    for (lmax, num, order) in vars:
        bspline_set = BSplineSet(order, 
                                 lin_knots(0.0, rmax, num))
        y_list = get_coupledY_set(lmax, 0, True, True)
        try:
            (h, s) = mat_h2(bond_length, bspline_set, y_list)
        except:
            print "Failed at building matrix"
            print "--------------------------------"
            print traceback.format_exc(sys.exc_info()[2])
            print "--------------------------------"
        write_coo_mat(h, "hmat.dat")
        write_coo_mat(s, "smat.dat")
    t1 = time.clock()
    ru = resource.getrusage(resource.RUSAGE_SELF)
    with open("mat.dat", "w") as f:
        f.write("time: {0}\n".format(t1-t0))
        f.write("mem(Mb): {0}\n".format(ru.ru_maxrss/(1024.0)))
    os.chdir("..")
    
