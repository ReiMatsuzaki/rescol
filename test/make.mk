test_angmoment.out:test_angmoment.o angmoment.o y1s.o y2s.o
test_bps.out: test_bps.o bps.o
test_bspline.out: test_bspline.o bspline.o bps.o pot.o mat.o cscaling.o viewerfunc.o
test_dvr.out: test_dvr.o dvr.o bps.o mat.o viewerfunc.o
#test_fd.out: test_fd.o fd.o bps.o mat.o
test_fem_inf.out: test_fem_inf.o fem_inf.o dvr.o bspline.o fd.o bps.o pot.o cscaling.o mat.o viewerfunc.o
test_mat.out: test_mat.o mat.o viewerfunc.o
test_oce2.out: test_oce2.o oce2.o mat.o fem_inf.o dvr.o bspline.o fd.o bps.o pot.o cscaling.o mat.o angmoment.o viewerfunc.o
test_pot.out: test_pot.o pot.o viewerfunc.o
test_cscaling.out: test_cscaling.o cscaling.o
TEST_LIST=angmoment bps bspline dvr fem_inf mat oce2 pot cscaling #fd

check_%: test_%.out
	@echo
	@echo ////Unit Test////
	./$<
	@echo
