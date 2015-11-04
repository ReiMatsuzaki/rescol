include ../common/common.mk
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

TEST_BIN_DIR=${BIN_DIR}/test
TEST_NAMES=angmoment bps bspline dvr fem_inf mat oce2 pot scale
TESTS=$(addprefix utest_, ${TEST_NAMES})
LIBGTEST=-lgtest

${OBJ_DIR}/%.o : ./%.c
	${PCC} -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} -o $@ $< 

${OBJ_DIR}/%.o : ./%.cpp
	${CXX} -c ${PCC_FLAGS} ${CCPPFLAGS} -o $@ $< 

utest_%: ${OBJ_DIR}/test_%.o ${LIB_OBJS}	
	@echo 
	@echo //////$@////////
	@-${CLINKER} -o ${TEST_BIN_DIR}/$@.out $^ \
	${LIBGSL} ${LIBGTEST} ${SLEPC_EPS_LIB}
	${TEST_BIN_DIR}/$@.out

all_test: ${TESTS}
	@echo ${TESTS}
