#include <slepceps.h>
#include <gtest/gtest.h>
extern "C" {
#include "../src_c/oce2.h"
}

static char help[] = "Unit test for oce_two.c";

TEST(TestFirst, First) {
  ASSERT_EQ(2, 1+1);
}

class TestOCE2 : public ::testing::Test {
public:
  BPS bps;
  BSS bss;
  FEMInf fem;
  Y2s y2s;
  OCE2 oce2;

  virtual void SetUp() {
    MPI_Comm comm = MPI_COMM_SELF;
    
    BPSCreate(&bps, comm); BPSSetExp(bps, 30.0, 21, 5.0);
    BSSCreate(&bss, 2, bps, comm);
    FEMInfCreateBSS(&fem, bss);
    Y2sCreate(&y2s, comm); Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);

    OCE2Create(&oce2, MPI_COMM_SELF); 
    OCE2Set(oce2, fem, y2s); 
  }
  virtual void TearDown() {
    OCE2Destroy(&oce2);
  }

};

TEST_F(TestOCE2, He) {

  Mat H, S;
  
  OCE2SetTMat(oce2, &H);
  OCE2PlusVneMat(oce2, 0.0, 1.0, &H);
  OCE2PlusVeeMat(oce2, &H);
  OCE2SetSMat(oce2, &S);

  EPS eps;
  EPSCreate(oce2->comm, &eps);  
  EPSSetProblemType(eps, EPS_GHEP);  
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 
  EPSSetTarget(eps, -4.0);  
  EPSSetOperators(eps, H, S);
  EPSSetFromOptions(eps); 

  EPSSolve(eps);

  PetscInt n;
  EPSGetConverged(eps, &n);
  ASSERT_TRUE(n > 0);

  PetscScalar kr;
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  
  ASSERT_NEAR(kr, -2.858256, 0.000001);
}
TEST_F(TestOCE2, view) {
  OCE2View(oce2);
}

int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}
