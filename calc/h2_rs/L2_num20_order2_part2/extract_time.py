from var import *

with open("res_time.dat", "w") as f_out:
    def func(bond_length):
        res_mat = keyval_to_dict("mat.out.dat")
        res_diag = keyval_to_dict("diag_detail.dat")
        f_out.write("{0} {1} {2} {3} {4}\n"
                    .format(bond_length, 
                            res_mat["t_calc"], 
                            res_mat["t_write"], 
                            res_diag["t_read"],
                            res_diag["t_diag"]
                    ))
        
    cd ("calc.out")
    for (dir_name, bond_length) in bond_length_list:
        with_dir(dir_name, lambda : func(bond_length))

