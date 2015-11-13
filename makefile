include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
# include ${RESCOL_DIR}/common/variable.mk
RESCOL_DIR ?= $(realpath .)

INC_DIR=${RESCOL_DIR}/include

VPATH = ${RESCOL_DIR}/src/lib \
	${RESCOL_DIR}/src/bin \
	${RESCOL_DIR}/test \
	${RESCOL_DIR}/test/training

CPPFLAGS+=-I${INC_DIR} -I${PETSC_DIR}
CFLAGS+=-std=gnu99
CXXFLAGS=
LIBS=-lgsl -lgtest

OBJ_FEM= fem_inf.o fd.o bspline.o dvr.o bps.o mat.o pot.o cscaling.o viewerfunc.o

include $(subst .c,.d,${SOURCE_LIST})
include ${RESCOL_DIR}/test/make.mk
include ${RESCOL_DIR}/test/training/make.mk
include ${RESCOL_DIR}/src/bin/make.mk


${BIN_DIR}:
	mkdir -p ${BIN_DIR}
${TEST_BIN_DIR}:
	mkdir -p ${TEST_BIN_DIR}

%.o: %.c ${OBJ_DIR}
	@echo [compile] $(notdir $<)
	@${PCC} -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} ${CPPFLAGS} -o $@ $< 
%.o : %.cpp ${OBJ_DIR}
	@echo [compile] $(notdir $<)
	@${CXX} -c ${PCC_FLAGS}  ${CXXFLAGS} ${CCPPFLAGS} ${CPPFLAGS} -o $@ $< 
%.d: %.c
	${PCC} -M ${CPPFLAGS} ${CCPPFLAGS}$< > $@.$$$$;                  \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.d: %.cpp
	${CXX} -M ${CPPFLAGS} ${CCPPFLAGS}$< > $@.$$$$;                  \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.out : %.o
	@echo [  link  ] $^
	@${CLINKER} -o $@ $^ ${LIBS} ${SLEPC_EPS_LIB}

.PHONY: TAGS
TAGS: 
	find ${PETSC_DIR}/include -regex '.*\.h$$' | grep -v ftn-auto | gtags -q -f -
	find ${PETSC_DIR}/${PETSC_ARCH}/include -regex '.*\.h$$' | grep -v ftn-auto | gtags -q -f -

.PHONY: check
check: $(addprefix check_, ${TEST_LIST})

.PHONY:myrm
myrm: 
	rm -f ${OBJ_DIR}/*.o
	rm -f ${BIN_DIR}/*.out
	rm -f ${TEST_BIN_DIR}/*.out
	rm -f ./*.o
	rm -f ./*.out

