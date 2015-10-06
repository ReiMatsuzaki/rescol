import sys
import os
import commands
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time
import resource
import traceback

bond_length = 1.4

t0 = time.clock()
try:
    os.system("../../../src_c/diag.out -eps_nev 1 -eps_real -eps_target -3.0")
except:
    print "Failed at diagonalizing matrix"
    print "--------------------------------"
    print traceback.format_exc(sys.exc_info()[2])
    print "--------------------------------"
    
t1 = time.clock()
print "time: ", t1-t0
ru = resource.getrusage(resource.RUSAGE_SELF)
print "mem(Mb): ", ru.ru_maxrss/(1024.0)

