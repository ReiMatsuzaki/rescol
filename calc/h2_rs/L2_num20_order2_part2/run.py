import os

def cd(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    os.chdir(dir_name)

script_dir = os.environ["RESCOL_DIR"] + "/script"
orig_dir = os.path.abspath(".")
bond_length_list = [ (str(i+3), 0.1*(i+3)) for i in range(4)]

cd("calc")
for (dir_name, bond_length) in bond_length_list:

    cd(dir_name)
    os.system("sed s/xxx/{0}/ {1}/mat.in > mat.in.dat".format(bond_length,
                                                                  orig_dir))
    os.system("cp {0}/diag.in .".format(orig_dir))
    os.system("python {0}/calc_mat.py".format(script_dir))
    #    os.system("{0}/calc_diag.py".format(script_dir))
    cd("..")

