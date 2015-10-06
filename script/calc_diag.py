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
os.system("{0}/src_c/diag.out".format(os.environ["RESCOL_DIR"]))
t1 = time.clock()

with open("diag.out.dat", "w") as f:
    f.write("t_calc: {0}\n".format(t1-t0))


