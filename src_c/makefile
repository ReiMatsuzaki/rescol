#$CXXFLAGS=-g
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

test_hatom: hatom.o bspline.o
	-${CLINKER} -o hatom.out  hatom.o bspline.o ${SLEPC_EPS_LIB}
	./hatom.out -bss_order 6 -bss_rmax 30.0 -bss_knots_num 100 -eps_nev 5

test_mat: mat.o test_mat.o chkopts
	-${CLINKER} -o test_mat.out  test_mat.o mat.o ${SLEPC_EPS_LIB}
	./test_mat.out -h-mat hmat.dat

diag.out: diag.o mat.o chkopts
	-${CLINKER} -I/usr/local/Cellar/petsc/3.6.1_3/real/include/petsc/mpiuni/ -o $@ diag.o mat.o ${SLEPC_EPS_LIB}

test_diag: diag.out chkopts
	python test_diag.py

h2plus.out: h2plus.o chkopts
	-${CLINKER} -o $@ h2plus.o bspline.o mat.o ${SLEPC_EPS_LIB}

test_h2plus: h2plus.out chkopts
	./h2plus.out -target_dir h2plus -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 2.0 -eps_target -1.0 -eps_nev 1 -eps_max_it 1000
	echo Reference: −1.102 634 214 494 9

test_h2plus_h_L0: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L0 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -4.0 -eps_nev 10 -eps_max_it 1000

test_h2plus_h_L1: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L1 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -0.51 -eps_nev 10 -eps_max_it 1000

test_h2plus_h: h2plus.out chkopts
	./h2plus.out -target_dir h2plus -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -2.2 -eps_nev 10 -eps_max_it 1000

h2mole.out: h2mole.o chkopts
	-${CLINKER} -o $@ h2mole.o bspline.o mat.o ${SLEPC_EPS_LIB}

test_h2mole: h2mole.out chkopts
	python ../script/calc_y2mat.py -t h2mole -l1 2
	./h2mole.out -target_dir h2mole -bss_order 6 -bss_knots_num 21 -bss_rmax 10.0 -qmax 10 -bond_length 0.0 -eps_target -2.2 -eps_nev 1 -eps_max_it 1000

