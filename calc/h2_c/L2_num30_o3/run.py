import os

param = {"order": 4, "num": 41, "rmax": 80.0}

x_list = range(24)
for x in x_list:
    os.system("rm -r {0}".format(x))


os.system("rm -r common")
os.system("mkdir common")
os.chdir("common")
os.system("""
python ../../../../script/calc_y2mat.py \
-l0 0 \
-l1 2 \
-l1guess 0 \
-l2guess 0 | tee calc_y2mat.out
""")
os.chdir("..")

os.system("rm -r 0")
os.system("mkdir -p 0")
script = """
../../../src_c/he_guess.out -L1 0 -L2 0 \
			    -z 1.44 \
			    -in_dir common \
			    -out_dir 0 \
			    -bss_knots_type exp \
			    -bss_rmax {rmax} \
			    -bss_order {order}\
			    -bss_knots_num {num} | tee he_guess.out.dat
""".format(**param)
sh_path = "0/run_he_guess.sh"
with open(sh_path, "w") as f:
    f.write(script)
os.system("sh "+sh_path)

for x in x_list:
    print x
    os.system("mkdir -p {0}".format(x))
    os.system("cp common/* {0}".format(x))
    if(x != x_list[0]):
        os.system("cp {0}/eig0.vec.dat {1}/".format(x-1, x))
        os.system("cp {0}/eig0.vec.dat {0}/guess.vec.dat".format(x))
    
    x_param = {"dir": x, "length": x/10.0}
    x_param.update(param)
    script = """
    ../../../src_c/h2mole2.out -in_dir {dir} \
    -out_dir {dir} \
    -bss_order {order} \
    -bss_knots_num {num} \
    -bss_knots_type exp \
    -bss_rmax {rmax} \
    -eri direct \
    -qmax 10 \
    -bond_length {length} \
    -guess_type vec \
    -eps_type jd \
    -eps_target -4.2 \
    -eps_nev 1 | tee {dir}/h2mole.out.dat \
    """.format(**x_param)
    path = "{0}/run_h2mole.sh".format(x)
    with open(path, "w") as f:
        f.write(script)
    os.system("sh " + path)



