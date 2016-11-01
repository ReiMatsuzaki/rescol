../../../bin/driv_1d \
    -L 1  \
    -Z 1.0 \
    -S_pot "sto 2.0 2 1.0" \
    -fem_type dvr \
    -bps_num_zs 201 -bps_zmax 100.0 -bps_type line \
    -dvr_nq 5 \
    -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 40.0 \
    -energy_range 0.5 \
    -viewerfunc_path psi.csv \
    -viewerfunc_num 501 -viewerfunc_xmax 0.25 -viewerfunc_xmax 100.25

echo ""
echo "reference alpha = (5.65688402161-1.08811622008j)"
