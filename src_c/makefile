CFLAGS=-std=gnu99
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

test_mat: mat.o test_mat.o
	-${CLINKER} -o test_mat.out  test_mat.o mat.o ${SLEPC_EPS_LIB}
	./test_mat.out -h-mat hmat.dat

hatom.out: hatom.o bspline.o mat.o chkopts

diag.out: diag.o mat.o chkopts
	-${CLINKER} -I/usr/local/Cellar/petsc/3.6.1_3/real/include/petsc/mpiuni/ -o $@ diag.o mat.o ${SLEPC_EPS_LIB}

h2plus.out: h2plus.o bspline.o mat.o chkopts
	-${CLINKER} -o $@ h2plus.o bspline.o mat.o ${SLEPC_EPS_LIB}

h2mole.out: h2mole.o bspline.o mat.o chkopts
	-${CLINKER} -o $@ h2mole.o bspline.o mat.o ${SLEPC_EPS_LIB}

h2mole2.out: h2mole2.o bspline.o mat.o chkopts
	-${CLINKER} -o $@ h2mole2.o bspline.o mat.o ${SLEPC_EPS_LIB}

he_guess.out: he_guess.o mat.o bspline.o chkopts
	-${CLINKER} -o $@ he_guess.o bspline.o mat.o ${SLEPC_EPS_LIB}

test_hatom: hatom.o bspline.o mat.o
	-${CLINKER} -o hatom.out  hatom.o bspline.o mat.o ${SLEPC_EPS_LIB}
	./hatom.out \
	-bss_order 6 \
	-bss_rmax 30.0 \
	-bss_knots_num 61 \
	-bss_knots_type exp \
	-eps_nev 5

test_diag: diag.out chkopts
	python test_diag.py

test_h2plus: h2plus.out chkopts
	./h2plus.out \
	-target_dir h2plus \
	-bss_order 4 \
	-bss_knots_num 31 \
	-bss_knots_type line \
	-bss_rmax 30.0 \
	-qmax 32 \
	-bond_length 2.0 \
	-eps_target -1.0 \
	-eps_nev 1 \
	-eps_max_it 1000
	echo Reference: 1.102 634 214 494 9

test_h2plus_h_L0: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L0 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -4.0 -eps_nev 10 -eps_max_it 1000

test_h2plus_h_L1: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L1 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -0.51 -eps_nev 10 -eps_max_it 1000

test_h2plus_h: h2plus.out chkopts
	./h2plus.out -target_dir h2plus -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -2.2 -eps_nev 10 -eps_max_it 1000

# this test takes 145s+82s=4min
test_h2mole: h2mole.out chkopts
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./h2mole.out -in_dir tmp -out_dir tmp -bss_order 4 -bss_knots_num 31 -bss_rmax 20.0 -bss_knots_type exp -qmax 10 -bond_length 0.0 -eps_target -3.0 -eps_nev 1 -eps_max_it 1000

# this test takes 
# 15s+82s = 2min (default eigen solver)
# 15s+197s = 4min (lanczos)
# 15s+157s =      (arnoldi)
# 15s+2s   (gd, -2.90375)
test_h2mole2: h2mole2.out
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./h2mole2.out -in_dir tmp -out_dir tmp -bond_length 0.0 -eps_target -3.0  -bss_order 4 -bss_knots_num 31 -bss_rmax 30.0 -bss_knots_type exp -qmax 10 -eps_max_it 1000 -eps_target -5.0 -eps_nev 1 -eps_type krylovschur
