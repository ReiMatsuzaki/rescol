# see CWMcCurdy and MFernando, JPB 37, (2004), 917
# check convergence speed

bin=../../../../../build/complex/driv1d.out
for r0 in 50 60 70 80 90
do
    echo $r0
    $bin -fem_type bss -bss_order 8 \
	 -bps_num_zs 101 -bps_zmax 100.0 -bps_type line \
	 -cscaling_type sharp_ecs -cscaling_r0 $r0 -cscaling_theta 20.0 \
	 -energy_range 0.5 \
	 -pot_type single \
	 -problem_type scatter \
	 -L 0 \
	 -v0-pot "sto -3.0 0 1.0" \
	 -viewerfunc_view ascii:psi_${r0}.dat \
	 -viewerfunc_num 2000 -viewerfunc_xmax 100.0 > res_${r0}.dat
done

