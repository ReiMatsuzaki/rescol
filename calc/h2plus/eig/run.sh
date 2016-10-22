../../../bin/eig_h2plus \
    -R 2.0 \
    -fem_type dvr \
    -dvr_nq 6 \
    -bps_num_zs 81 -bps_zmax 40.0 -bps_type line \
    -y1s_lmax 6 \
    -eps_target -3.5
#    -viewerfunc_view ascii:psi.dat \
#    -viewerfunc_num 501 -viewerfunc_xmax 100.0

echo "reference energy is -1.1026342144949."
echo "from H.Wind, JCP 42, 2371 (1965)"
#python -c "print -16.2*3.674931*0.01"
#python -c "print -16.2*3.674931*0.01+1.0/2.0"


    
