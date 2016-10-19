include local.mk
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

OBJDIR=$(realpath obj)
BINDIR=$(realpath bin)

VPATH=lib:test

OBJS=$(SRCS_C:.c=.o) $(SRCS_CPP:.cpp=.o)
DEPS=$(SRCS_C:.c=.d) $(SRCS_CPP:.cpp=.d)

# -- google test --
# read README in googletest
GTEST_DIR=$(HOME)/local/src/googletest/googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $(OBJDIR)/$@
gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $(OBJDIR)/$@ $(OBJDIR)/$^

.SUFFIXES:
.SUFFIXES: .o .c .cpp

%.o : %.cpp
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $(OBJDIR)/$@
%.o : %.c
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $(OBJDIR)/$@

test_bps.out : test_bps.o gtest.a bps.o
	cd $(OBJDIR); ${CXX} -o $(BINDIR)/$@ $^  ${SLEPC_EPS_LIB}

check_%: test_%.out
	$(BINDIR)/$<

.PHONY: clean
clean::
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	rm -f $(BINDIR)/*.out

