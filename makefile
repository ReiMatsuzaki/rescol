include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
include local.mk

## verbosity controll
### see http://stackoverflow.com/questions/9314322/controlling-verbosity-of-make
V = 0
AT_0:=@
AT_1:=
AT=$(AT_$(V))
PCC_:=$(PCC)
PCC=$(AT)$(PCC_)
CXX_:=$(CXX)
CXX=$(AT)$(CXX_)

OBJDIR=./obj
BINDIR=./bin
VPATH=src:lib:test
SRCS_c:=$(wildcard src/*.c) $(wildcard lib/*.c) $(wildcard test/*.cpp) \
SRCS_cpp:=$(wildcard test/*.cpp) \
DEPS:=$(SRCS_cpp:%.cpp=${OBJDIR}/%.d) $(SRCS_c:%.c=${OBJDIR}/%.d)

-include $(DEPS)

# -- google test --
# read README in googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

$(OBJDIR)/gtest-all.o : $(GTEST_SRCS_)
	@echo "build" $(notdir $@)
	$(AT)$(CC) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@
$(OBJDIR)/gtest.a : $(OBJDIR)/gtest-all.o
	@echo "build" $(notdir $@)
	$(AT)$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/%.o : %.cpp
	@echo "build" $(notdir $@)
	$(AT)[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(AT)[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	$(CXX) $(PCC_FLAGS) $(CCPPFLAGS) $(CXXFLAGS) $(CPPFLAGS) -MMD -c $< -o $@

$(OBJDIR)/%.o : %.c
	@echo "build" $(notdir $@)
	$(AT)[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(AT)[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	$(PCC) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -MMD -c $< -o $@


OBJS = test_bps.o bps.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_bps.out : $(OBJSFULL)
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB}

OBJS = test_angmoment.o y1s.o y2s.o angmoment.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_angmoment.out : $(OBJSFULL)
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_cscaling.o cscaling.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_cscaling.out : $(OBJSFULL)
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_mat.o mat.o synthesize.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_mat.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_pot.o pot.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_pot.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_synthesize.o synthesize.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_synthesize.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_bspline.o bspline.o cscaling.o bps.o eeps.o pot.o mat.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_bspline.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_dvr.o dvr.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_dvr.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_fem_inf.o fd.o fem_inf.o dvr.o bspline.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o viewerfunc.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_fem_inf.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_oce1.o oce1.o fem_inf.o fd.o dvr.o bspline.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o viewerfunc.o y1s.o angmoment.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_oce1.out : $(OBJSFULL)
	@echo "build" $@
	${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

check_%: $(BINDIR)/test_%.out
	@echo $@
	$(AT)$<

.PHONY: check
check:: check_bps check_angmoment check_cscaling check_mat check_pot check_synthesize check_bspline check_dvr check_fem_inf check_oce1

.PHONY: clean
clean::
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	rm -f $(BINDIR)/*.out

