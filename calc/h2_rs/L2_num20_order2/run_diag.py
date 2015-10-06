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

os.chdir("calc")
for bond_length in bond_length_list:
    dir_name = "{0}".format(round(bond_length*10))
    print dir_name
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    os.chdir(dir_name)
    t0 = time.clock()
    try:
        os.system("../../../../../src_c/diag.out -eps_nev 1 -eps_target_real -eps_target -3.0")
    except:
        print "Failed at diagonalizing matrix"
        print "--------------------------------"
        print traceback.format_exc(sys.exc_info()[2])
        print "--------------------------------"
    
    t1 = time.clock()
    ru = resource.getrusage(resource.RUSAGE_SELF)
    mem = ru.ru_maxrss/(1024.0)
    with open("diag.dat", "w") as f:
        f.write("time: {0}\n".format(t1-t0))
        f.write("mem: {0}\n")
    os.chdir("..")

