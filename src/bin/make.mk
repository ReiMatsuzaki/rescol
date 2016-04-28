# ==== plot basis ====
plot_basis.out: plot_basis.o ${OBJ_FEM}


# ==== Solve Driv Eq ====
.PHONY: check_driv1d_dvr
check_driv1d_dvr: driv1d.out
	@for num in 201 401 601; do \
		echo $$num; \
		./$< -fem_type dvr -dvr_nq 8 \
		-bps_num_zs $$num -bps_zmax 100.0 -bps_type line \
		-energy_range 0.5 \
		-driv-pot "sto -2.0 2 1.0" \
		-pot_type single \
		-problem_type driv \
		-L 1 \
		-cscaling_r0 70.0 -cscaling_theta 20.0 \
		-v0-pot "pow -1.0 -1" \
		-viewerfunc_view ascii:tmp/driv1d_dvr.dat -viewerfunc_num 100 -viewerfunc_xmax 100.0 | grep c0;\
	done
	@echo "ref(a): 5.65688402161, 1.08811622008"

check_driv1d_dvr_r0: driv1d.out
	@for R0 in 50.0 100.0 15.0 170.0; do \
		./$< -fem_type dvr -dvr_nq 8 \
		-bps_num_zs 201 -bps_zmax 200.0 -bps_type line \
		-energy_range 0.5 \
		-driv-pot "sto -2.0 2 1.0" \
		-pot_type single \
		-problem_type driv \
		-L 1 \
		-cscaling_r0 $$R0 -cscaling_theta 20.0 \
		-v0-pot "pow -1.0 -1" \
		-viewerfunc_view ascii:tmp/driv1d_dvr.dat -viewerfunc_num 100 -viewerfunc_xmax 100.0 | grep c0;\
	done
	@echo "ref(a):    5.65688402161, -1.08811622008"

# see 2016/4/16
driv1d.out: driv1d.o ${OBJ_FEM}
.PHONY: check_driv1d
check_driv1d: driv1d.out
	@for num in 201 401 601; do \
		echo $$num; \
		./$< -fem_type bss -bss_order 8 \
		-bps_num_zs $$num -bps_zmax 200.0 -bps_type line \
		-cscaling_type sharp_ecs -cscaling_r0 150.0 -cscaling_theta 5.0 \
		-energy_range 0.5 \
		-L 1 \
		-driv-pot "sto -2.0 2 1.0" \
		-pot_type single \
		-problem_type driv \
		-v0-pot "pow -1.0 -1" \
		-viewerfunc_view ascii:tmp/driv1d.dat \
		-viewerfunc_num 100 -viewerfunc_xmax 100.0 | grep c0; \
	done
#Reference from 2016/4/12/* and test_r1gtoint
	@echo "ref(a): 5.65688402161, 1.08811622008"

.PHONY: check_driv1d2
check_driv1d2: driv1d.out
	./$< -fem_type bss -bss_order 8 \
	-bps_num_zs 501 -bps_zmax 100.0 -bps_type line \
	-cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
	-energy_range 0.3:0.5:2 \
	-driv-pot "sto -2.0 2 1.0" \
	-pot_type double \
	-problem_type driv \
	-L 1 \
	-v0-pot "pow -1.0 -1" \
	-v1-pot "sto -1.0 0 1.0" \
	-viewerfunc_view ascii:tmp/driv1d.dat -viewerfunc_num 100 -viewerfunc_xmax 100.0

check_driv1d_scatter: driv1d.out
	./$< -fem_type bss -bss_order 5 \
	-bps_num_zs 401 -bps_zmax 100.0 -bps_type line \
	-cscaling_type sharp_ecs -cscaling_r0 50.0 -cscaling_theta 5.0 \
	-energy_range 0.5 \
	-pot_type single \
	-problem_type scatter \
	-L 0 \
	-v0-pot "sto -3.0 0 1.0" \
	-viewerfunc_view ascii:tmp/driv1d.dat \
	-viewerfunc_num 100 -viewerfunc_xmax 100.0 | grep cross
	@echo "ref(cross sec) = 10.108998"
#@echo "ref(phase) = 2.028 859 730 528"

install_driv1d: driv1d.out
	cp driv1d.out ~/bin/driv1d


h_pi.out: h_pi.o ${OBJ_FEM} y1s.o angmoment.o

eps_help.out:

h2mole.out: h2mole.o ${OBJ_FEM} oce2.o angmoment.o eeps.o y2s.o synthesize.o y1s.o
basis=-fem_type bss -fem_type bss -bss_order 2 -bps_num_zs 21 -bps_zmax 20.0 -bps_type exp -y2s_rot sigma -y2s_parity gerade -y2s_mirror plus -y2s_lmax 0
.PHONY: check_h2mole_bss
check_h2mole_bss: h2mole.out he_guess.out
	./he_guess.out $(basis) \
	-guess_type fit -z 1.5 -L1 0 -L2 0 \
	-vec_viewer binary:guess.dat
	./h2mole.out $(basis) \
	-eps_nev 1 -eps_max_it 1000 -eps_type jd -eps_view_values ascii -eps_target -4.0 \
	-init_path guess.dat \
	-bondlength 0.0 -num_bondlength 1 -d_bondlength 0.5 -malloc_dump
#	./h2mole.out $(basis) \
	-eps_nev 2 -eps_max_it 1000 -eps_type jd -eps_view_values ascii -eps_target -4.0 \
	-init_path guess2.dat \
	-bondlength 0.0 -num_bondlength 1 -d_bondlength 0.5 -malloc_dump

he_guess.out: he_guess.o ${OBJ_FEM} oce2.o angmoment.o y2s.o y1s.o eeps.o synthesize.o wavefunc.o
.PHONY: check_he_guess
check_he_guess: he_guess.out 
	-rm tmp/*
	./he_guess.out \
	-fem_type bss -fem_type bss -bss_order 2 -bps_num_zs 21 -bps_zmax 20.0 -bps_type exp \
	-guess_type fit -z 1.0 -L1 0 -L2 0 \
	-vec_viewer binary:guess.dat -malloc_dump


.PHONY: check_h_pi
check_h_pi: h_pi.out
	./$< -fem_type bss -bss_order 4 \
	-bps_num_zs 101 -bps_zmax 100.0 -bps_type line \
	-cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
	-viewerfunc_view ascii:tmp/h_pi.dat -viewerfunc_num 100 -viewerfunc_xmax 100.0
	head tmp/h_pi.dat
install_h_pi: h_pi.out
	cp h_pi.out ~/bin/h_pi

fit.out: ${OBJ_FEM} viewerfunc.o pot.o synthesize.o wavefunc.o
.PHONY: check_fit
check_fit: fit.out
	./$< -fem_type bss -bss_order 4 \
	-bps_num_zs 51 -bps_zmax 51.0 -bps_type line \
	-wavefunc_type heig -wavefunc_heig_n 1 -wavefunc_heig_l 0 \
	-viewerfunc_xmax 20.0 -viewerfunc_num 200 -viewerfunc_view ascii:tmp/fit.dat


write_pot.out: ${OBJ_FEM} oce1.o angmoment.o pot.o viewerfunc.o
.PHONY: check_write_pot
check_write_pot: write_pot.out
	./$< -pot_type slater -pot_v0 3.5 -pot_z 1.0 \
	-viewerfunc_view ascii:stdout -viewerfunc_num 20 -viewerfunc_xmax 20.0 \
	-malloc_dump


eig_one.out: eig_one.o ${OBJ_FEM} oce1.o angmoment.o y1s.o eeps.o viewerfunc.o synthesize.o
.PHONY: check_eig_one
check_eig_one: eig_one.out
	./$< -fem_type bss -bss_order 4 \
	-bps_num_zs 51 -bps_zmax 51.0 -bps_type line \
	-cscaling_type sharp_ecs -cscaling_r0 40.0 -cscaling_theta 25.0 \
	-y1s_L 0 \
	-pot_type slater -pot_v0 3.5 -pot_n 2 -pot_z 1.0 \
	-eps_nev 2 \
	-eps_converged_reason \
	-eeps_view_values ::ascii_info_detail \
	-viewerfunc_xmax 100.0 -viewerfunc_num 10 -viewerfunc_view ascii:stdout \
	-malloc_dump

.PHONY: check_eig_one_h
check_eig_one_h: eig_one.out
	./$< -fem_type bss -bss_order 4 -bps_num_zs 21 -bps_zmax 20.0 -bps_type exp \
	-y1s_L 0 \
	-pot_type coulomb -pot_z -1.0 \
	-eps_nev 2 -eps_target -1.0 -eps_view_values ascii \
	-malloc_dump


h2plus.out: h2plus.o ${OBJ_FEM} oce1.o y1s.o angmoment.o eeps.o
.PHONY: check_h2plus
check_h2plus: h2plus.out
	./$< -fem_type bss -bss_order 8 \
	-bps_num_zs 101 -bps_zmax 50.0 -bps_type line \
	-y1s_rot sigma -y1s_parity gerade -y1s_lmax 8 \
	-eps_nev 2 -eps_max_it 1000 -eps_target -1.2 -eeps_view_values ascii:stdout
	@echo Reference: -1.102 634 214 494 9

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

