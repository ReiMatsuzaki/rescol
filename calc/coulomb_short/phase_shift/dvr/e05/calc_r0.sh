for r0 in 70 80 90 100 110 120 130 140 150 160 170 180 190 200
do
    echo $(( $r0 +31))
    driv1d -fem_type dvr -dvr_nq 8 \
	   -bps_num_zs $(( $r0 +31)) -bps_zmax $(($r0+30)) -bps_type line \
	   -cscaling_type sharp_ecs -cscaling_r0 ${r0} -cscaling_theta 20.0 \
	   -pot_type double \
	   -problem_type driv \
	   -energy_range 0.5 \
	   -driv-pot "sto -2.0 2 1.0" \
	   -L 1 \
	   -Z 1.0 \
	   -v1-pot "sto -1.0 0 1.0" \
	   > res_r0_${r0}.dat    
done

sh get_r0.sh
python plot_r0.py
