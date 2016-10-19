include local.mk
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

OBJDIR=$(realpath obj)
BINDIR=$(realpath bin)

VPATH=src:lib:test

SRCS_c:=$(wildcard src/*.c) $(wildcard lib/*.c) $(wildcard test/*.cpp) \
SRCS_cpp:=$(wildcard test/*.cpp) \
DEPS:=$(SRCS_cpp:%.cpp=${DIR}/%.d) $(SRCS_c:%.c=${DIR}/%.d)
-include ${DEPS}

# -- google test --
# read README in googletest
GTEST_DIR=$(HOME)/local/src/googletest/googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

$(OBJDIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@
$(OBJDIR)/gtest.a : $(OBJDIR)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/%.o : %.cpp
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MMD -c $< -o $@
$(OBJDIR)/%.o : %.c
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -MMD -c $< -o $@

OBJS = test_bps.o bps.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_bps.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB}

OBJS = test_angmoment.o y1s.o y2s.o angmoment.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_angmoment.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_cscaling.o cscaling.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_cscaling.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_mat.o mat.o synthesize.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_mat.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_pot.o pot.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_pot.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_synthesize.o synthesize.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_synthesize.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_bspline.o bspline.o cscaling.o bps.o eeps.o pot.o mat.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_bspline.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_dvr.o dvr.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_dvr.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_fem_inf.o fd.o fem_inf.o dvr.o bspline.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o viewerfunc.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_fem_inf.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

OBJS = test_oce1.o oce1.o fem_inf.o fd.o dvr.o bspline.o cscaling.o bps.o eeps.o pot.o synthesize.o mat.o viewerfunc.o y1s.o angmoment.o gtest.a
OBJSFULL=$(foreach o, $(OBJS), $(OBJDIR)/$o)
$(BINDIR)/test_oce1.out : $(OBJSFULL)
	cd $(OBJDIR); ${CXX} -o $@ $^  ${SLEPC_EPS_LIB} -lgsl

check_%: $(BINDIR)/test_%.out
	$<

.PHONY: clean
clean::
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	rm -f $(BINDIR)/*.out

