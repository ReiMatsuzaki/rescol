TESTS=angmoment bps bspline dvr fem_inf mat oce2 pot scale
GO_TESTS=$(addprefix go_test_, ${TESTS})

${TEST_BIN_DIR}/test_bspline.out: ${OBJ_DIR}/test_bspline.o ${OBJ_DIR}/bspline.o ${OBJ_DIR}/bps.o ${OBJ_DIR}/pot.o ${OBJ_DIR}/mat.o ${OBJ_DIR}/scale.o
	@echo ////link test_bspline////
	@-${CLINKER} -o $@ $^ ${LIBGSL} ${LIBGTEST} ${SLEPC_EPS_LIB}

${TEST_BIN_DIR}/test_%.out: ${OBJ_DIR}/test_%.o ${LIB_OBJS}	
	@echo link $(notdir $<)
	@-${CLINKER} -o $@ $^ ${LIBGSL} ${LIBGTEST} ${SLEPC_EPS_LIB}

#.PHONY: go_test_bspline
#go_test_bspline: ${TEST_BIN_DIR}/test_bspline.out
#	@echo 
#	@echo //////$@////////
#	@-./$<

.PHONY: go_test_%
go_test_%: ${TEST_BIN_DIR}/test_%.out
	@echo 
	@echo //////$@////////
	@-./$<

.PHONY: all_test
all_test: ${GO_TESTS}
