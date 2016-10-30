#!/bin/zsh
bindir=../../../bin
fem=(-fem_type dvr -dvr_nq 6 -bps_num_zs 201 -bps_zmax 100.0 -bps_type line)
lmax=7
target=(-R 2.0)

${bindir}/fit_oce1 ${fem} -y1s_lmax ${lmax} \
	 -v_pot "sto 2.0 1 1.0" \
	 -out fit.dat\
	 -info fit_oce1_info -log_trace fit_oce1_trace -malloc_dump

${bindir}/eig_h2plus ${target} ${fem} -y1s_lmax ${lmax} \
    -eps_target -3.5 \
    -in fit.dat \
    -out cs0.dat \
    -viewerfunc_view ascii:psi0.dat \
    -viewerfunc_num 501 -viewerfunc_xmax 100.0 \
    -info eig_h2plus_info -log_trace eig_h2plus_trace -malloc_dump

${bindir}/driv_h2plus ${target} ${fem} \
	 -init_y1s_lmax  ${lmax} -init_y1s_rot  sigma -init_y1s_parity gerade\
	 -final_y1s_lmax ${lmax} -final_y1s_rot sigma -final_y1s_parity ungerade\
	 -cscaling_type sharp_ecs -cscaling_r0 70.0 -cscaling_theta 40.0 \
	 -viewerfunc_path psi1.dat -viewerfunc_num 1000 -viewerfunc_xmax 100.0 \
	 -in cs0.dat \
	 -out cs1.dat \
	 -E0 -1.099966 \
	 -w_range 2.099966 \
	 ## reference is 0.0165619128315(GTO) or 0.015962(ref) at E=1.0

