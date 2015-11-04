FEMOBJS=../src_c/fem_inf.o ../src_c/fd.o ../src_c/bspline.o ../src_c/dvr.o ../src_c/bps.o ../src_c/mat.o ../src_c/pot.o ../src_c/scale.o
ANGOBJS=../src_c/angmoment.o

CFLAGS=-std=gnu99
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

test_mat:test_mat.o
	-${CLINKER} -o $@.out $@.o ../src_c/mat.o ${SLEPC_EPS_LIB}
	./$@.out

test_angmoment:test_angmoment.o
	-${CLINKER} -o test_angmoment.out test_angmoment.o ../src_c/angmoment.o -lgtest -lgsl ${SLEPC_EPS_LIB}
	./test_angmoment.out

test_bspline:test_bspline.o ../src_c/bspline.o
	-${CLINKER} -o test_bspline.out test_bspline.o ${FEMOBJS} ${SLEPC_EPS_LIB} -g
	./test_bspline.out

test_dvr:test_dvr.o ../src_c/dvr.o ../src_c/bps.o
	-${CLINKER} -o test_dvr.out test_dvr.o ../src_c/dvr.o ../src_c/bps.o ../src_c/mat.o -llapack -lblas -lm ${SLEPC_EPS_LIB} -g
	./test_dvr.out -eps_max_it 1000 -eps_type jd

test_bps: test_bps.o ../src_c/bps.o 
	-${CLINKER} -o test_bps.out test_bps.o ../src_c/bps.o ${SLEPC_EPS_LIB} -g
	./test_bps.out -bps_type exp -bps_num_zs 6 -bps_zmax 2.0

test_fem_inf: test_fem_inf.o ../src_c/fem_inf.o
	-${CLINKER} -o test_fem_inf.out test_fem_inf.o ${FEMOBJS} ${SLEPC_EPS_LIB} -g
	./test_fem_inf.out

test_oce2: test_oce2.o ../src_c/fem_inf.o ../src_c/oce2.o
	-${CLINKER} -o test_oce2.out test_oce2.o ../src_c/oce2.o ${FEMOBJS} ${ANGOBJS} -lgtest -lgsl ${SLEPC_EPS_LIB} -g
	./test_oce2.out -petsc_dump

test_scale.out: test_scale.o ../src_c/scale.o
	-${CLINKER} -o $@ test_scale.o ../src_c/scale.o -lgtest ${SLEPC_EPS_LIB} -g

test_scale: test_scale.out
	./test_scale.out -petsc_dump
