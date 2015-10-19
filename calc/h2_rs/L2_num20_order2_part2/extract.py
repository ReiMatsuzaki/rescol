from var import *

def extract_eigvals(bond_length, f_out):
    with open("diag_eigvals.dat", "r") as f_in:
        f_in.readline() # skip header
        line_list = f_in.readlines()
        ene_list = [float(line.split(" ")[0]) for line in line_list] 
        f_out.write("{0} ".format(bond_length))
        for ene in ene_list:
            f_out.write("{0} ".format(ene+1.0/bond_length))
        f_out.write("\n")

with open("res.dat", "w") as f_out:
    cd("calc.out")
    for (dir_name, bond_length) in bond_length_list:
        with_dir(dir_name, 
                 lambda: extract_eigvals(bond_length, f_out)
                 )
