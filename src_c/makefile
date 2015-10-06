#$CXXFLAGS=-g
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

test_mat: mat.o test_mat.o chkopts
	-${CLINKER} -o test_mat.out  test_mat.o mat.o ${SLEPC_EPS_LIB}
	./test_mat.out -h-mat hmat.dat

diag.out: diag.o mat.o chkopts
	-${CLINKER} -o $@ diag.o mat.o ${SLEPC_EPS_LIB}

test_diag: diag.out chkopts
	python test_diag.py
