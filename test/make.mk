test_angmoment.out:test_angmoment.o angmoment.o y1s.o y2s.o
test_bps.out: test_bps.o bps.o
test_bspline.out: test_bspline.o bspline.o bps.o pot.o mat.o cscaling.o viewerfunc.o eeps.o
test_dvr.out: test_dvr.o dvr.o bps.o mat.o viewerfunc.o synthesize.o
#test_fd.out: test_fd.o fd.o bps.o mat.o
test_fem_inf.out: test_fem_inf.o fem_inf.o dvr.o bspline.o fd.o bps.o pot.o cscaling.o mat.o viewerfunc.o eeps.o synthesize.o
test_mat.out: test_mat.o mat.o viewerfunc.o
test_oce2.out: test_oce2.o oce2.o mat.o fem_inf.o dvr.o bspline.o fd.o bps.o pot.o cscaling.o mat.o angmoment.o viewerfunc.o synthesize.o y2s.o y1s.o
test_oce1.out: test_oce1.o oce1.o mat.o fem_inf.o dvr.o bspline.o fd.o bps.o pot.o cscaling.o mat.o angmoment.o y1s.o viewerfunc.o eeps.o synthesize.o
test_pot.out: test_pot.o pot.o viewerfunc.o
test_cscaling.out: test_cscaling.o cscaling.o
test_op.out: test_op.o
test_synthesize.out: test_synthesize.o synthesize.o
TEST_LIST=angmoment bps bspline dvr fem_inf mat oce2 pot cscaling synthesize
#fd op

default_opt=-fem_type bss -bss_order 4 -bps_num_zs 40 -bps_zmax 100.0 -bps_type line
time_synthesize.out: time_synthesize.o ${OBJ_FEM} synthesize.o
time_synthesize: time_synthesize.out
	./time_synthesize.out ${default_opt} -mat_mat_synthesize_alg 0
	./time_synthesize.out ${default_opt} -mat_mat_synthesize_alg 2
	./time_synthesize.out ${default_opt} -mat_mat_synthesize_alg 3

check_%: test_%.out
	@echo
	@echo ////Unit Test////
	./$< -malloc_dump
	@echo
