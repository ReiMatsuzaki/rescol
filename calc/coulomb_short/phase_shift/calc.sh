for order in 8
do
    for num in 500
    do
	echo $order $num
	driv1d -fem_type bss -bss_order $order \
	       -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	       -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
	       -energy_range 0.1:1.0:100 \
	       -pot_type double \
	       -driv-pot "sto -2.0 2 1.0" \
	       -L 1 \
	       -problem_type driv \
	       -v0-pot "pow -1.0 -1|pow 1.0 -2" \
	       -v1-pot "sto -1.0 0 1.0" \
	       -viewerfunc_view ascii:short.dat -viewerfunc_num 500 -viewerfunc_xmax 100.0 > short_res.dat   > res_${order}_${num}.dat    
    done
done
