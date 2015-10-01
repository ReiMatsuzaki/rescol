# petsc and slepc settings
SLEPC_DIR = ${HOME}/local/src/slepc-3.4.3
PETSC_DIR = ${HOME}/local/src/petsc-3.4.3
PETSC_ARCH = rccls_C
#$CXXFLAGS=-g

test_mat.out: mat.o test_mat.o chkopts
	-${CLINKER} -o $@  test_mat.o mat.o ${SLEPC_LIB}
	./test_mat.out -h-mat hmat.dat

diag.out: diag.o mat.o chkopts
	-${CLINKER} -o $@ diag.o mat.o ${SLEPC_LIB}	

include ${SLEPC_DIR}/conf/slepc_common
