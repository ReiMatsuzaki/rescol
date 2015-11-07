h2mole.out: h2mole.o ${OBJ_FEM} oce2.o angmoment.o
he_guess.out: he_guess.o ${OBJ_FEM} oce2.o angmoment.o
h_pi.out: h_pi.o ${OBJ_FEM} angmoment.o

.PHONY: check_h2mole_bss
check_h2mole_bss: he_guess.out h2mole.out
	rm -f tmp/*
	./he_guess.out -z 1.5 -in_dir tmp -out_dir tmp -guess_type calc \
	-fem_type bss -bss_order 2 -bps_num_zs 21 -bps_type exp -bps_zmax 30.0
	@echo 
	./h2mole.out -eri direct -guess_type read \
	-in_dir tmp -out_dir tmp  \
	-fem_type bss -bss_order 2 -bps_num_zs 21 -bps_type exp -bps_zmax 30.0 \
	-y2s_rot sigma -y2s_parity gerade -y2s_mirror plus -y2s_lmax 2 \
	-eps_nev 1 -eps_max_it 1000 -eps_type jd \
	-bondlength 0.0 -num_bondlength 1 -d_bondlength 0.5

.PHONY: check_he_guess
check_he_guess: he_guess.out
	-rm tmp/*
	./he_guess.out -z 1.0 \
	-in_dir tmp -out_dir tmp \
	-guess_type eig \
	-fem_type fd -fd_num 100 -fd_xmax 20.0

.PHONY: check_h_pi
check_h_pi: h_pi.out
	./$< -fem_type bss -bss_order 4 \
	-bps_num_zs 101 -bps_zmax 100.0 -bps_type line \
	-scaler_type secs -scaler_r0 70.0 -scaler_theta 20.0



check_hatom: hatom.o fem_inf.o fd.o dvr.o bspline.o mat.o
	-${CLINKER} -o hatom.out  hatom.o fem_inf.o fd.o dvr.o bspline.o bps.o mat.o ${SLEPC_EPS_LIB}
	./hatom.out -fem_type fd -fd_num 1000 -fd_xmax 100.0 -eps_nev 5	 \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 -guess_type heig
	./hatom.out -fem_type bss -bss_order 6 -eps_nev 5 \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 
	./hatom.out -fem_type dvr -dvr_nq 5 -eps_nev 5 -eps_type jd \
	-bps_num_zs 31 -bps_type exp -bps_zmax 30.0 

check_diag: diag.out chkopts
	python test_diag.py

run_h2plus: h2plus.out chkopts
	./h2plus.out \
	-y1s_rot sigma -y1s_parity gerade -y1s_lmax 14 \
	-fem_type bss -bss_order 8 -bps_num_zs 41 -bps_type exp -bps_zmax 30.0 \
	-bond_length 2.0 -eps_nev 1 -eps_max_it 1000 -eps_type jd
	echo Reference: 1.102 634 214 494 9

run_h2plus_h_L0: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L0 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -4.0 -eps_nev 10 -eps_max_it 1000

test_h2plus_h_L1: h2plus.out chkopts
	./h2plus.out -target_dir h2plus_L1 -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -0.51 -eps_nev 10 -eps_max_it 1000

test_h2plus_h: h2plus.out chkopts
	./h2plus.out -target_dir h2plus -bss_order 6 -bss_knots_num 61 -bss_rmax 30.0 -qmax 32 -bond_length 0.0 -eps_target -2.2 -eps_nev 10 -eps_max_it 1000



test_h2mole_fd: h2mole.out chkopts
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./he_guess.out -z 1.5 -in_dir tmp -out_dir tmp -guess_type calc \
	-fem_type fd -fd_xmax 30.0 -fd_num 300 
	./h2mole.out -bond_length 0.0 -eri direct -guess_type read  -qmax 32 \
	-in_dir tmp -out_dir tmp  \
	-fem_type fd -fd_xmax 30.0 -fd_num 300 \
	-eps_nev 1 -eps_max_it 1000 -eps_type jd \
	-bondlength0 0.0 -num_bondlength 2 -d_bondlength 0.5

test_h2mole_old_bss: h2mole.out chkopts
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./he_guess.out -z 1.5 -in_dir tmp -out_dir tmp -guess_type calc \
	-fem_type bss -bss_order 4 -bps_num_zs 21 -bps_type exp -bps_zmax 30.0
	./h2mole.out -bond_length 0.0 -eri direct -guess_type read  -qmax 32 \
	-in_dir tmp -out_dir tmp  \
	-fem_type bss -bss_order 4 -bps_num_zs 21 -bps_type exp -bps_zmax 30.0 \
	-eps_nev 1 -eps_max_it 1000 -eps_type jd \
	-bondlength0 0.0 -num_bondlength 2 -d_bondlength 0.5

test_h2mole_dvr: h2mole.out chkopts
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./he_guess.out -z 1.0 -in_dir tmp -out_dir tmp -guess_type calc \
	-fem_type dvr -dvr_nq 4 -bps_num_zs 11 -bps_type exp -bps_zmax 30.0
	./h2mole.out -bond_length 0.0 -eri direct -guess_type read  -qmax 32 \
	-in_dir tmp -out_dir tmp  \
	-fem_type dvr -dvr_nq 4 -bps_num_zs 11 -bps_type exp -bps_zmax 30.0 \
	-eps_nev 1 -eps_max_it 1000 -eps_type krylovschur \
	-bondlength0 0.0 -num_bondlength 2 -d_bondlength 0.5


# this test takes 
# 15s+82s = 2min (default eigen solver)
# 15s+197s = 4min (lanczos)
# 15s+157s =      (arnoldi)
# 15s+2s   (gd, -2.90375)
test_h2mole2: h2mole2.out
	rm tmp/*
	python ../script/calc_y2mat.py -t tmp -l1 2 > tmp/y2mat.dat
	./h2mole2.out -in_dir tmp -out_dir tmp -bond_length 0.0 -eps_target -3.0  -bss_order 4 -bss_knots_num 31 -bss_rmax 30.0 -bss_knots_type exp -qmax 10 -eps_max_it 1000 -eps_target -5.0 -eps_nev 1 -eps_type krylovschur

