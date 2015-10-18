import sys
import time
from os.path import abspath, dirname, join
sys.path.append(join(dirname(dirname(__file__)), "src"))
from rescol import *

try:
    key_val = keyval_to_dict("diag.in.dat")
except:
    print "failed to open file: diag.in.dat"
    sys.exit()

nev = int(key_val["nev"])
target= float(key_val["target"])

t0 = time.clock()
os.system("{0}/src_c/diag.out -eps_nev {1} -eps_target_real -eps_target {2}"
          .format(os.environ["RESCOL_DIR"], nev, target))
t1 = time.clock()

