#/usr/bin/python

import time
import traceback
from os.path import abspath, dirname, join
import sys
sys.path.append(join(dirname(dirname(__file__)), "src"))
from rescol import *

"""
reading mat.in.dat and 
compute Hamiltonian and overlap matrix
write hmat.dat and smat.dat
"""

try:
    key_val = keyval_to_dict("mat.in.dat")
except:
    print "failed to open file: mat.in.dat"
    sys.exit()

bond_length = float(key_val["bond_length"])
order = int(key_val["order"])
rmax = float(key_val["rmax"])
num = int(key_val["num"])
lmax = int(key_val["lmax"])

bspline_set = BSplineSet(order, 
                         lin_knots(0.0, rmax, num))
y_list = get_coupledY_set(lmax, 0, True, True)

try:
    t0 = time.clock()
    (h, s) = mat_h2(bond_length, bspline_set, y_list)
    t1 = time.clock()
    t_calc = t1-t0
except:
    print "Failed at building matrix"
    print "--------------------------------"
    print traceback.format_exc(sys.exc_info()[2])
    print "--------------------------------"

t0 = time.clock()
write_coo_mat(h, "hmat.dat")
write_coo_mat(s, "smat.dat")
t1 = time.clock()
t_write = t1-t0

with open("mat.out.dat", "w") as f:
    f.write("t_calc: {0}\n".format(t_calc))
    f.write("t_write: {0}\n".format(t_write))
    f.write("n_basis: {0}\n".format(h.shape[0]))
    f.write("n_data: {0}\n".format(len(h.data)))


    
