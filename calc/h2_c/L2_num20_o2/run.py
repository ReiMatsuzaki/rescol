import os

x_list = range(25)
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
			    -bss_rmax 30.0 \
			    -bss_order 2\
			    -bss_knots_num 21 | tee he_guess.out.dat
"""
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
        
    script = """
    ../../../src_c/h2mole2.out -in_dir {0} \
			  -out_dir {0} \
			  -bss_order 2 \
			  -bss_knots_num 21 \
			  -bss_knots_type exp \
			  -bss_rmax 30.0 \
			  -eri direct \
			  -qmax 10 \
			  -bond_length {1} \
			  -guess_type vec \
			  -eps_target -4.2 \
			  -eps_nev 1 | tee {0}/h2mole.out.dat \
    """.format(x, x/10.0)
    path = "{0}/run_h2mole.sh".format(x)
    with open(path, "w") as f:
        f.write(script)
    os.system("sh " + path)



