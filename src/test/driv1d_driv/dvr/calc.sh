# see CWMcCurdy and MFernando, JPB 37, (2004), 917
# check convergence speed

bin=../../../../../build/complex/driv1d.out
for num in 101 201 301 401
do
    echo $num
    $bin -fem_type dvr -dvr_nq 8 \
	 -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	 -R0 70.0 -theta 20.0 \
	 -energy_range 0.5 \
	 -pot_type single \
	 -problem_type driv \
	 -L 1 \
	 -driv-pot "sto -2.0 2 1.0" \
	 -v0-pot "pow -1.0 -1" \
	 -viewerfunc_view ascii:psi_${num}.dat \
	 -viewerfunc_num 2000 -viewerfunc_xmax 100.0 > res_${num}.dat
done

