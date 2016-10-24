#!/bin/zsh
basis=(-fem_type dvr -dvr_nq 6 -bps_num_zs 81 -bps_zmax 40.0 -bps_type line -y1s_lmax 6)

echo "stat fit"
../../../bin/fit_oce1 ${basis} \
		      -v_pot "sto 2.0 1 1.0" \
		      -out fit.dat\
		      -info fit_oce1_info -log_trace fit_oce1_trace -malloc_dump

echo "stat eig"
../../../bin/eig_h2plus \
    -R 2.0 \
    ${basis} \
    -eps_target -3.5 \
    -in fit.dat \
    -out cs.dat \
    -viewerfunc_view ascii:psi.dat \
    -viewerfunc_num 501 -viewerfunc_xmax 100.0 \
    -log_summary -info eig_h2plus_info -log_trace eig_h2plus_trace -malloc_dump

echo "reference energy is -1.1026342144949."
echo "from H.Wind, JCP 42, 2371 (1965)"
#python -c "print -16.2*3.674931*0.01"
#python -c "print -16.2*3.674931*0.01+1.0/2.0"


    
