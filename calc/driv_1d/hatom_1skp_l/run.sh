../../../bin/driv_1d \
    -L 1  \
    -Z 1.0 \
    -S-pot "sto 2.0 2 1.0" \
    -fem_type dvr \
    -bps_num_zs 201 -bps_zmax 100.0 -bps_type line \
    -dvr_nq 5 \
    -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 20.0 \
    -energy-range 0.5 \
    -viewerfunc_view ascii:psi.dat \
    -viewerfunc_num 501 -viewerfunc_xmax 100.0

echo ""
echo "reference alpha = (5.65688402161-1.08811622008j)"
