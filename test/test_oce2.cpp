#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/oce2.h>

static char help[] = "Unit test for oce2.c";

class TestOCE2 : public ::testing::Test {
public:
  OCE2 oce2;

  virtual void SetUp() {
    MPI_Comm comm = MPI_COMM_SELF;
    
    BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 30.0, 21, 5.0);
    BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, 2, bps); BSSSetUp(bss);
    FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

    Y2s y2s; Y2sCreate(comm, &y2s); Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);

    OCE2Create(comm, &oce2); OCE2Set(oce2, fem, y2s);
    
  }
  virtual void TearDown() {
    OCE2Destroy(&oce2);
  }

};
TEST_F(TestOCE2, He) {

  PetscErrorCode ierr;
  if(getenv("SHOW_DEBUG")) {
    ierr = OCE2View(oce2, PETSC_VIEWER_STDOUT_SELF); ASSERT_EQ(0, ierr);
  }
  

  Mat H; 
  ierr = OCE2TMat(oce2, MAT_INITIAL_MATRIX, &H); ASSERT_EQ(0, ierr);
  ierr = OCE2PlusVneMat(oce2, 0.0, 1.0, H); ASSERT_EQ(0, ierr);
  ierr = OCE2PlusVeeMat(oce2, H); ASSERT_EQ(0, ierr);

  Mat S;  
  ierr = OCE2SMat(oce2, MAT_INITIAL_MATRIX, &S, NULL); ASSERT_EQ(0, ierr);

  EPS eps;
  EPSCreate(oce2->comm, &eps);  
  EPSSetProblemType(eps, EPS_GHEP);  
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 
  EPSSetType(eps, EPSJD);
  EPSSetTarget(eps, -4.2);
  EPSSetOperators(eps, H, S);
  EPSSetFromOptions(eps); 

  EPSSolve(eps);

  PetscInt n;
  EPSGetConverged(eps, &n);
  ASSERT_TRUE(n > 0);

  PetscScalar kr;
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  
  ASSERT_NEAR(PetscRealPart(kr), -2.858256, 0.000001);
}

/*
TEST(TestOCE2Dvr, Dvr) {
  
  MPI_Comm comm = MPI_COMM_SELF;
  
  BPS bps; BPSCreate(&bps, comm); BPSSetExp(bps, 30.0, 11, 5.0);
  DVR dvr; DVRCreate(&dvr, 3, bps, comm);
  FEMInf fem; FEMInfCreateDVR(&fem, dvr);
  Y2s y2s; Y2sCreate(&y2s, comm); Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);

  OCE2 oce2; OCE2Create(&oce2, MPI_COMM_SELF); 
  OCE2Set(oce2, fem, y2s); 

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
*/
int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}
