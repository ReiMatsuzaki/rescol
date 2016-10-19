# see CWMcCurdy and MFernando, JPB 37, (2004), 917
# check convergence speed

bin=../../../../../build/complex/driv1d.out
for num in 101 201 301 401
do
    $bin -fem_type dvr -dvr_nq 8 \
	 -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	 -R0 80.0 -theta 20.0 \
	 -energy_range 0.5 \
	 -pot_type single \
	 -problem_type scatter \
	 -L 0 \
	 -v0-pot "sto -3.0 0 1.0" \
	 -viewerfunc_view ascii:psi_${num}.dat \
	 -viewerfunc_num 2000 -viewerfunc_xmax 100.0 > res_${num}.dat
done

