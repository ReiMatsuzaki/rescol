include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

test_mat:test_mat.o
	-${CLINKER} -o $@.out $@.o ../src_c/mat.o ${SLEPC_EPS_LIB}
	./$@.out

test_angmoment:test_angmoment.o
	-${CLINKER} -o test_angmoment.out test_angmoment.o ../src_c/angmoment.o -lpython ${SLEPC_EPS_LIB}
	./test_angmoment.out

test_bspline:test_bspline.o ../src_c/bspline.o
	-${CLINKER} -o test_bspline.out test_bspline.o ../src_c/bspline.o ${SLEPC_EPS_LIB} -g
	./test_bspline.out







