rm *.dat
python ../../../script/calc_y2mat.py -l0 0 -l1 4 -l1guess 0 -l2guess 0 | tee calc_y2mat.out

../../../src_c/he_guess.out -L1 0 -L2 0 \
			    -z 1.0 \
			    -bss_knots_type exp \
			    -bss_rmax 30.0 \
			    -bss_order 4\
			    -bss_knots_num 21 | tee he_guess.out

../../../src_c/h2mole.out -bss_order 4 \
			  -bss_knots_num 21 \
			  -bss_knots_type exp \
			  -bss_rmax 30.0 \
			  -eri direct \
			  -qmax 10 \
			  -bond_length 0.0 \
			  -guess_type non_ee \
			  -eps_target -4.2 \
			  -eps_nev 1 \


