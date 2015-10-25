CFLAGS=-std=gnu99
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
FEMOBJS = fem_inf.o fd.o bspline.o dvr.o bps.o mat.o

test_mat: mat.o test_mat.o
	-${CLINKER} -o test_mat.out  test_mat.o mat.o ${SLEPC_EPS_LIB}
	./test_mat.out -h-mat hmat.dat

hatom.out: hatom.o fem_inf.o bspline.o dvr.o mat.o chkopts

diag.out: diag.o mat.o chkopts
	-${CLINKER} -I/usr/local/Cellar/petsc/3.6.1_3/real/include/petsc/mpiuni/ -o $@ diag.o mat.o ${SLEPC_EPS_LIB}

h2plus.out: h2plus.o ${FEMOBJS} chkopts
	-${CLINKER} -o $@ h2plus.o ${FEMOBJS} ${SLEPC_EPS_LIB}

h2mole.out: h2mole.o bspline.o mat.o chkopts
	-${CLINKER} -o $@ h2mole.o bspline.o mat.o ${SLEPC_EPS_LIB}

h2mole2.out: h2mole2.o fem_inf.o dvr.o bspline.o fd.o bps.o mat.o chkopts
	-${CLINKER} -o $@ h2mole2.o fem_inf.o dvr.o bspline.o fd.o bps.o mat.o ${SLEPC_EPS_LIB}

he_guess.out: he_guess.o mat.o bspline.o chkopts
	-${CLINKER} -o $@ he_guess.o bspline.o mat.o ${SLEPC_EPS_LIB}

he_guess_dvr.out: he_guess_dvr.o mat.o dvr.o chkopts
	-${CLINKER} -o $@ he_guess_dvr.o dvr.o mat.o ${SLEPC_EPS_LIB}

test_hatom: hatom.o fem_inf.o fd.o dvr.o bspline.o mat.o
	-${CLINKER} -o hatom.out  hatom.o fem_inf.o fd.o dvr.o bspline.o bps.o mat.o ${SLEPC_EPS_LIB}
	./hatom.out -fem_type fd -fd_num 1000 -fd_xmax 100.0 -eps_nev 5	 \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 -guess_type heig
	./hatom.out -fem_type bss -bss_order 6 -eps_nev 5 \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 
	./hatom.out -fem_type dvr -dvr_nq 5 -eps_nev 5 -eps_type jd \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 

test_diag: diag.out chkopts
	python test_diag.py

test_h2plus: h2plus.out chkopts
	./h2plus.out -target_dir h2plus -qmax 32  \
	-fem_type bss -bss_order 4 -bps_num_zs 21 -bps_type exp -bps_zmax 30.0 \
	-bond_length 2.0 -eps_nev 1 -eps_max_it 1000 -eps_type jd
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

test_h2mole_dvr: he_guess_dvr.out
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./he_guess_dvr.out -L1 0 -L2 0 -z 1.0 -in_dir tmp -out_dir tmp  -dvr_nq 5 -dvr_knots_num 10 -dvr_rmax 30.0 -dvr_knots_type exp -eps_type jd
