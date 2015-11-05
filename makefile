include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
RESCOL_DIR = .

SRC_DIR=${RESCOL_DIR}/src
LIB_SRC_DIR=${SRC_DIR}/lib
OBJ_DIR=${RESCOL_DIR}/obj/${PETSC_ARCH}
BIN_DIR=${RESCOL_DIR}/bin/${PETSC_ARCH}
INC_DIR=${RESCOL_DIR}/include

TEST_BIN_DIR=${BIN_DIR}/test

LIB_SRCS=fd.c scale.c bps.c fem_inf.c angmoment.c oce1.c oce2.c bspline.c dvr.c mat.c 	pot.c
LIB_OBJS=${addprefix ${OBJ_DIR}/, ${patsubst %.c,%.o,${LIB_SRCS}}}

VPATH = ${RESCOL_DIR}/src/lib ${RESCOL_DIR}/src/bin ${RESCOL_DIR}/test

CFLAGS=-std=gnu99 -I${INC_DIR}
CXXFLAGS= -I${INC_DIR}

LIBGSL=-lgsl
LIBGTEST = -lgtest

include src/bin/makefile
include test/makefile

${OBJ_DIR}/%.o : %.c ${OBJ_DIR}
	@echo compile $(notdir $<)
	${PCC} -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} -o $@ $< 
${OBJ_DIR}/%.o : %.cpp ${OBJ_DIR}
	@echo compile $(notdir $<)
	${CXX} -c ${PCC_FLAGS}  ${CXXFLAGS} ${CCPPFLAGS} -o $@ $< 

myrm: 
	rm -f ${OBJ_DIR}/*.o
	rm -f ${BIN_DIR}/*.out
	rm -f ${TEST_BIN_DIR}/*.out

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}
${BIN_DIR}:
	mkdir -p ${BIN_DIR}
${TEST_BIN_DIR}:
	mkdir -p ${TEST_BIN_DIR}

