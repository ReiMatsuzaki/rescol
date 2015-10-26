import os

rmax = 40.0
order = 4
num = 41
l1 = 4

x_list = range(25)
for x in x_list:
    os.system("rm -r {0}".format(x))

os.system("rm -r common")
os.system("mkdir common")
os.chdir("common")
os.system("""
python ../../../../script/calc_y2mat.py \
-l0 0 \
-l1 {0} \
-l1guess 0 \
-l2guess 0 | tee calc_y2mat.out
""".format(l1))
os.chdir("..")

os.system("rm -r 0")
os.system("mkdir -p 0")
os.system("""
../../../src_c/he_guess.out -L1 0 -L2 0 \
			    -z 1.00 \
			    -in_dir common \
			    -out_dir 0 \
			    -bss_knots_type exp \
			    -bss_rmax {0} \
			    -bss_order {1}\
			    -bss_knots_num {2} | tee he_guess.out.dat
""".format(rmax, order, num))

for x in x_list:
    print x
    os.system("mkdir -p {0}".format(x))
    os.system("cp common/* {0}".format(x))
    if(x != x_list[0]):
        os.system("cp {0}/eig0.vec.dat {1}/".format(x-1, x))
        os.system("cp {0}/eig0.vec.dat {0}/guess.vec.dat".format(x))
        
    os.system("""
    ../../../src_c/h2mole.out -in_dir {3} \
			  -out_dir {3} \
			  -bss_order {1} \
			  -bss_knots_num {2} \
			  -bss_knots_type exp \
			  -bss_rmax {0} \
			  -eri direct \
			  -qmax 10 \
			  -bond_length {4} \
			  -guess_type vec \
			  -eps_target -4.2 \
			  -eps_nev 1 | tee {3}/h2mole.out.dat \
    """.format(rmax, order, num, x, x/10.0))


