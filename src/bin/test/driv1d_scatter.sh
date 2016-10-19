# see CWMcCurdy and MFernando, JPB 37, (2004), 917
# check convergence speed

bin=../../../build/complex/driv1d.out
for num in 101 201 301 401
do
    $bin -fem_type bss -bss_order 8 \
	 -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	 -cscaling_type sharp_ecs -cscaling_r0 50.0 -cscaling_theta 40.0 \
	 -energy_range 0.5 \
	 -pot_type single \
	 -problem_type scatter \
	 -L 0 \
	 -v0-pot "sto -3.0 0 1.0" \
	 -viewerfunc_view ascii:tmp/driv1d.dat \
	 -viewerfunc_num 100 -viewerfunc_xmax 100.0 > res_${num}.dat
done
