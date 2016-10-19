bin=../../../../build/complex/plot_basis.out


for n in `seq 40 70`
do
    echo $n
    $bin -fem_type bss \
	 -bss_order 8 \
	 -bps_num_zs 71 \
	 -bps_zmax 70.0 \
	 -bps_type line \
	 -cscaling_type sharp_ecs \
	 -cscaling_r0 50.0 \
	 -cscaling_theta 40.0 \
	 -n ${n} \
	 -viewerfunc_view ascii:basis${n}.dat \
	 -viewerfunc_num 1000 -viewerfunc_xmax 100.0
done
