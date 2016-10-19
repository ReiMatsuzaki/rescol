for order in 4 6 8
do
    for num in 101 201 301
    do
	echo $order $num
	driv1d -fem_type dvr -dvr_nq $order \
	       -bps_num_zs $num -bps_zmax 100.0 -bps_type line \
	       -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
	       -pot_type double \
	       -problem_type driv \
	       -energy_range 0.5 \
	       -driv-pot "sto -2.0 2 1.0" \
	       -L 1 \
	       -Z 1.0 \
	       -v1-pot "sto -1.0 0 1.0" \
	       > res_${order}_${num}.dat    
    done
done
