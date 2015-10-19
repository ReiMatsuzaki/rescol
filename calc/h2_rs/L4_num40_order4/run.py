import os
import sys
from var import *

cd("calc.out")
for (dir_name, bond_length) in bond_length_list:
    with_dir(dir_name, 
             lambda : 
             os.system("sed s/xxx/{0}/ {1}/mat.in > mat.in.dat"
                       .format(bond_length, orig_dir))
             or os.system("cp {0}/diag.in diag.in.dat".format(orig_dir))
             or os.system("python {0}/calc_mat.py".format(script_dir))
             or os.system("python {0}/calc_diag.py".format(script_dir))    
             )
