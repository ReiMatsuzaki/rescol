#include <gtest/gtest.h>
#include <slepceps.h>
#include "../include/oce1.h"
#include "../include/eeps.h"

static char help[] = "Unit test for oce1.c";

class TestOCE1 :public ::testing::Test {
public:
  MPI_Comm comm;
  OCE1 oce;   
  virtual void SetUp() {
    comm = PETSC_COMM_SELF;
    Y1s y1s; Y1sCreate(comm, &y1s); Y1sSetOne(y1s, 0, 1);
    BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
    int order = 5;
    BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
    FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);
    OCE1Create(comm, &this->oce); OCE1Set(this->oce, fem, y1s);    
  }
  virtual void TearDown() {
    OCE1Destroy(&this->oce);
  }
};
TEST_F(TestOCE1, HAtom) {

  if(getenv("SHOW_DEBUG"))
    OCE1View(oce, PETSC_VIEWER_STDOUT_SELF);

  Pot pot; PotCreate(comm, &pot); PotSetCoulombNE(pot, 0, 0.0, -1.0);
  Mat S; OCE1SMat(oce, MAT_INITIAL_MATRIX, &S, NULL);
  Mat H; OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  Mat V; OCE1PotMat(oce, ROT_SCALAR, pot, MAT_INITIAL_MATRIX, &V);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);
  
  EEPS eps; EEPSCreate(comm, &eps); 
  EEPSSetOperators(eps, H, S);
  EEPSSetTarget(eps, -0.2);
  EEPSSolve(eps);

  PetscScalar k;
  EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);
  ASSERT_NEAR(-0.125, PetscRealPart(k), 0.00001);

  PFDestroy(&pot);
  MatDestroy(&S); MatDestroy(&H); MatDestroy(&V);
  EEPSDestroy(&eps);
}
TEST_F(TestOCE1, HAtom2) {
  Pot pot; PotCreate(comm, &pot); PotSetCoulombNE(pot, 0, 0.0, -1.0);

  Mat S; OCE1SMat(oce, MAT_INITIAL_MATRIX, &S, NULL);
  Mat H; OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  OCE1PlusPotMat(oce, ROT_SCALAR, pot, H);
  
  EEPS eps; EEPSCreate(comm, &eps); 
  EEPSSetOperators(eps, H, S);
  EEPSSetTarget(eps, -0.2);
  EEPSSolve(eps);

  PetscScalar k;
  EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);
  ASSERT_NEAR(-0.125, PetscRealPart(k), 0.00001);

  PFDestroy(&pot);
  MatDestroy(&S); MatDestroy(&H);
  EEPSDestroy(&eps);
  
}

TEST(TestOCE1DVR, HAtom) {
  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s y1s;
  Y1sCreate(comm, &y1s); Y1sSet(y1s, SIGMA, GERADE, 4);
  BPS bps;
  BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
  DVR dvr;
  DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps); DVRSetUp(dvr);
  FEMInf fem;
  FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);
  OCE1 oce;
  OCE1Create(comm, &oce); OCE1Set(oce, fem, y1s);
  
  Pot pot;
  PotCreate(comm, &pot); PotSetCoulombNE(pot, 0, 0.0, -1.0);
  Mat H; OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  Mat V; OCE1PotMat(oce, ROT_SCALAR, pot, MAT_INITIAL_MATRIX, &V);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);
  
  EEPS eps;
  ierr = EEPSCreate(comm, &eps); ASSERT_EQ(0, ierr);
  ierr = EEPSSetOperators(eps, H, NULL); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps, -0.2); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps); ASSERT_EQ(0, ierr);

  PetscScalar k;
  EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);
  ASSERT_NEAR(-0.125, PetscRealPart(k), 0.0001);

  PFDestroy(&pot);
  MatDestroy(&H); MatDestroy(&V);
  EEPSDestroy(&eps);
}
TEST(TestOce1DVR, H2plus) {

  PetscErrorCode ierr;

  PetscReal R = 2.0;

  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s y1s;
  Y1sCreate(comm, &y1s); Y1sSet(y1s, SIGMA, GERADE, 8);
  BPS bps;
  BPSCreate(comm, &bps); BPSSetLine(bps, 40.0, 41);
  DVR dvr;
  DVRCreate(comm, &dvr); DVRSetKnots(dvr, 6, bps); DVRSetUp(dvr);
  FEMInf fem;
  FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);
  OCE1 oce;
  OCE1Create(comm, &oce); OCE1Set(oce, fem, y1s);
  
  Mat H; OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  OCE1PlusVneMat(oce, R/2.0, 1.0, H);

  // initial guess
  KSP ksp;
  ierr = KSPCreate(comm, &ksp); ASSERT_EQ(0, ierr);
  Pot slater;
  ierr = PotCreate(comm, &slater); ASSERT_EQ(0, ierr);
  ierr = PotSetSlater(slater, 3.0, 1, 2.0); ASSERT_EQ(0, ierr);
  Vec c_guess;
  ierr = OCE1CreateVec(oce, &c_guess); ASSERT_EQ(0, ierr);
  ierr = OCE1Fit(oce, slater, 0, ksp, c_guess); ASSERT_EQ(0, ierr);  
  
  EEPS eps;
  ierr = EEPSCreate(comm, &eps); ASSERT_EQ(0, ierr);
  ierr = EPSSetInitialSpace(eps->eps, 1, &c_guess); ASSERT_EQ(0, ierr);  
  ierr = EEPSSetOperators(eps, H, NULL); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps, -3.0); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps); ASSERT_EQ(0, ierr);

  PetscScalar k;
  ierr = EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);
  ASSERT_EQ(0, ierr);
  ASSERT_NEAR(-1.1026342144949, PetscRealPart(k), 0.002);

  OCE1Destroy(&oce); MatDestroy(&H);
  KSPDestroy(&ksp)   PFDestroy(&slater); 
  VecDestroy(&c_guess); 
  EEPSDestroy(&eps);  

}
/*
class TestOCE1H2plus :public ::testing::Test {
public:
  MPI_Comm comm;
  OCE1 oce;   
  virtual void SetUp() {
    comm = PETSC_COMM_SELF;

    Y1s y1s; Y1sCreate(comm, &y1s); Y1sSet(y1s, SIGMA, GERADE, 4);
    BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 40.0, 41);
    int order = 4;
    BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
    FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);
    OCE1Create(comm, &this->oce); OCE1Set(this->oce, fem, y1s);    
  }
  virtual void TearDown() {
    OCE1Destroy(&this->oce);
  }
};
TEST_F(TestOCE1H2plus, H2plusMat) {

  PetscErrorCode ierr;
  OceH2plus h2plus;
  ierr = OCE1CreateH2plus (oce, 0.3, 1.0, &h2plus); ASSERT_EQ(0, ierr);

  int nr, ny; OCE1GetSizes(oce, &nr, &ny);

  Mat H, S, H0, S0;
  ierr = OCE1H2plusMat(oce, h2plus, &H, &S, NULL);  ASSERT_EQ(0, ierr);
  ierr = OCE1H2plusMat_direct(oce, h2plus, &H0, &S0, NULL);  ASSERT_EQ(0, ierr);

  PetscScalar *val; PetscMalloc1(nr*ny, &val);
  for(int i = 0; i < nr*ny; i++)
    val[i] = i*0.2;
  Vec x;  VecCreateSeqWithArray(comm, nr*ny, nr*ny, val, &x);
  Vec y;  VecCreate(comm, &y); VecSetSizes(y, PETSC_DECIDE, nr*ny); VecSetUp(y);
  Vec y0; VecCreate(comm, &y0);VecSetSizes(y0, PETSC_DECIDE, nr*ny); VecSetUp(y0);

  MatMult(H, x, y);
  MatMult(H0, x, y0);

  VecAXPY(y, -1.0, y0);
  PetscReal norm;
  VecNorm(y, NORM_1, &norm);
  ASSERT_NEAR(0.0, norm/(nr*ny), pow(10.0,-10.0));
  
  VecDestroy(&x); VecDestroy(&y); VecDestroy(&y0);
  MatDestroy(&H); MatDestroy(&S); MatDestroy(&H0); MatDestroy(&S0);
  ierr = OCE1H2plusDestroy(&h2plus);  ASSERT_EQ(0, ierr);
  PetscFree(val);
}
TEST_F(TestOCE1H2plus, H2plusEPS) {

  PetscErrorCode ierr;
  OceH2plus h2plus;
  ierr = OCE1CreateH2plus (oce, 1.0, 1.0, &h2plus); ASSERT_EQ(0, ierr);
  Mat H, S;
  ierr = OCE1H2plusMat(oce, h2plus, &H, &S, NULL);  ASSERT_EQ(0, ierr);

  EEPS eps; 
  ierr = EEPSCreate(comm, &eps);  ASSERT_EQ(0, ierr);
  ierr = EEPSSetOperators(eps, H, S); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps, -3.2); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps); ASSERT_EQ(0, ierr);

  PetscScalar k;
  ierr = EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);  ASSERT_EQ(0, ierr);
  ASSERT_NEAR(-1.1026342144949, PetscRealPart(k), 0.03); 
  
  ierr = EEPSDestroy(&eps); ASSERT_EQ(0, ierr);
  ierr = OCE1H2plusDestroy(&h2plus);  ASSERT_EQ(0, ierr);
  ierr = MatDestroy(&H); ASSERT_EQ(0, ierr);
  ierr = MatDestroy(&S); ASSERT_EQ(0, ierr);
}
TEST_F(TestOCE1H2plus, H2plusEPS_direct) {

  PetscErrorCode ierr;
  OceH2plus h2plus;
  ierr = OCE1CreateH2plus (oce, 1.0, 1.0, &h2plus); ASSERT_EQ(0, ierr);
  Mat H, S;
  ierr = OCE1H2plusMat_direct(oce, h2plus, &H, &S, NULL);  ASSERT_EQ(0, ierr);

  EEPS eps; 
  ierr = EEPSCreate(comm, &eps);  ASSERT_EQ(0, ierr);
  ierr = EEPSSetOperators(eps, H, S); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps, -3.2); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps); ASSERT_EQ(0, ierr);

  PetscScalar k;
  ierr = EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);  ASSERT_EQ(0, ierr);
  ASSERT_NEAR(-1.1026342144949, PetscRealPart(k), 0.03); 
  
  ierr = EEPSDestroy(&eps); ASSERT_EQ(0, ierr);
  ierr = OCE1H2plusDestroy(&h2plus);  ASSERT_EQ(0, ierr);
  ierr = MatDestroy(&H); ASSERT_EQ(0, ierr);
  ierr = MatDestroy(&S); ASSERT_EQ(0, ierr);
}
*/
int _main(int argc, char **args) {
 ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  _main(argc, args);
  SlepcFinalize();
  return 0;
}
