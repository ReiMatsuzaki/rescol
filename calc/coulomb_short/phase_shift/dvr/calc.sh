for order in 8
do
    for num in 101
    do
	echo $order $num
	driv1d -fem_type dvr -dvr_nq $order \
	       -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	       -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
	       -energy_range 0.01:2.00:200 \
	       -pot_type double \
	       -problem_type driv\
	       -Z 1.0 \
	       -L 1 \
	       -driv-pot "sto -2.0 2 1.0" \
	       -v1-pot "sto -1.0 0 1.0" \
	       -viewerfunc_view ascii:short.dat -viewerfunc_num 500 -viewerfunc_xmax 100.0 > short_res.dat > res_${order}_${num}.dat    
    done
done

sh get.sh
