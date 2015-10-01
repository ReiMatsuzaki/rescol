import sys
import commands
sys.path.append("../../../src")
from rescol import *
import scipy.sparse.linalg as la
import time
import resource
import traceback

bond_length = 1.4
h = read_coo_mat("h.out")
s = read_coo_mat("s.out")

t0 = time.clock()
try:
    (es, vecs) = la.eigsh(h, 1, s, sigma=-3.0)
except:
    print "Failed at diagonalizing matrix"
    print "--------------------------------"
    print traceback.format_exc(sys.exc_info()[2])
    print "--------------------------------"
    
t1 = time.clock()
print "E0: ", es[0]+1.0/bond_length
print "time: ", t1-t0
ru = resource.getrusage(resource.RUSAGE_SELF)
print "mem(Mb): ", ru.ru_maxrss/(1024.0)

