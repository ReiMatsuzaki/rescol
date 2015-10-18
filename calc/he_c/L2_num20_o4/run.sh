rm *.dat
python ../../../script/calc_y2mat.py -l0 0 -l1 2 | tee calc_y2mat.out
../../../src_c/h2mole.out -bss_order 4 \
			  -eri direct \
			  -bss_knots_num 21 \
			  -bss_knots_type exp\
			  -bss_rmax 30.0 \
			  -qmax 10 \
			  -bond_length 0.0 \
			  -eps_target -4.2 \
			  -eps_nev 1 \
			  -eps_max_it 1000  | tee h2mole.out
