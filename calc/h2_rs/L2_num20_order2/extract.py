import os

bond_length_list = [i*0.1+0.3 for i in range(40)]

res = open("res.dat", "w")
os.chdir("calc")
for bond_length in bond_length_list:
    dir_name = "{0}".format(round(bond_length*10))
    print dir_name
    if not os.path.exists(dir_name):
        os.mkdir(dir_name) 
    os.chdir(dir_name)
    with open("diag_eigvals.dat", "r") as f:
        f.readline()
        [ene, err] = map(float, f.readline().split(" "))
        res.write("{0} {1} {2}\n".format(bond_length, ene, ene+1.0/bond_length))
    os.chdir("..")

res.close()
