import os
import sys

src_dir =os.environ["RESCOL_DIR"] + "/src"
script_dir = os.environ["RESCOL_DIR"] + "/script"
orig_dir = os.path.abspath(".")
bond_length_list = [ (str(i+3), 0.1*(i+3)) for i in range(4)]

sys.path.append(src_dir)
from rescol import cd, with_dir

cd("calc.out")
for (dir_name, bond_length) in bond_length_list:
    with_dir(dir_name, 
             lambda : 
             os.system("sed s/xxx/{0}/ {1}/mat.in > mat.in.dat".format(bond_length,
                                                                       orig_dir))
             or os.system("cp {0}/diag.in .".format(orig_dir))
             or os.system("python {0}/calc_mat.py".format(script_dir))
             # or os.system("{0}/calc_diag.py".format(script_dir))    
             )
