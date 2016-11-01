#include <gtest/gtest.h>
#include <slepceps.h>
#include "../include/oce1.h"
#include "../include/eeps.h"

static char help[] = "Unit test for oce1.c";

PetscErrorCode VecMatVecMult(Vec a, Mat M, Vec b, PetscScalar *res) {

  PetscErrorCode ierr;
  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)M, &comm); CHKERRQ(ierr);

  int na, nM1, nM2, nb;
  VecGetSize(a, &na);
  MatGetSize(M, &nM1, &nM2);
  VecGetSize(b, &nb);

  if(na != nM1 || nb != nM2) {
    SETERRQ(comm, 1, "size mismatch");
  }
  
  Vec Mb;
  ierr = MatCreateVecs(M, &Mb, NULL); CHKERRQ(ierr);
  ierr = MatMult(M, b, Mb); CHKERRQ(ierr);
  
  PetscScalar aMb;
  ierr = VecTDot(a, Mb, &aMb); CHKERRQ(ierr);
  *res = aMb;

  VecDestroy(&Mb);
  
  return 0;
}

class TestOCE1 :public ::testing::Test {
public:
  MPI_Comm comm;
  OCE1 oce;   
  virtual void SetUp() {
    comm = PETSC_COMM_SELF;
    Y1s y1s; Y1sCreate(comm, &y1s); Y1sSetOne(y1s, 1, 0);
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
  ierr = EEPSSetTarget(eps, -0.6); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps); ASSERT_EQ(0, ierr);

  PetscScalar k;
  Vec c; MatCreateVecs(H, NULL, &c);
  EEPSGetEigenpair(eps, 0, &k, c);
  ASSERT_NEAR(-0.5, PetscRealPart(k), 0.0001);

  PetscScalar norm;
  VecTDot(c, c, &norm);
  ASSERT_NEAR(1.0, PetscRealPart(norm), 0.0000001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(norm), 0.0000001);
  
  PetscScalar x = 2.5;
  PetscScalar y_ref = 2.0 * x * exp(-x);
  PetscScalar y_calc;
  OCE1PsiOne(oce, c, 0, 0, x, &y_calc);
  ASSERT_NEAR(PetscRealPart(y_ref),
	      PetscRealPart(y_calc), 0.0001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(y_calc), 0.0000001);

  OCE1Destroy(&oce); PFDestroy(&pot);
  MatDestroy(&H); MatDestroy(&V);
  VecDestroy(&c);
  EEPSDestroy(&eps);
}
TEST(TestOCE1DVR, HAtom_2p) {
  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s y1s;
  Y1sCreate(comm, &y1s); Y1sSet(y1s, 0, UNGERADE, 5);
  BPS bps;
  BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
  DVR dvr;
  DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps); DVRSetUp(dvr);
  FEMInf fem;
  FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);
  OCE1 oce;
  OCE1Create(comm, &oce); OCE1Set(oce, fem, y1s);
  // OCE1View(oce, PETSC_VIEWER_STDOUT_SELF);
  
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
  Vec c; MatCreateVecs(H, NULL, &c);
  EEPSGetEigenpair(eps, 0, &k, c);

  PetscScalar norm;
  VecTDot(c, c, &norm);
  ASSERT_NEAR(1.0, PetscRealPart(norm), 0.0000001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(norm), 0.0000001);

  PetscScalar x = 2.5;
  PetscScalar y_ref, y_calc;
  y_ref = 1.0/sqrt(24.0) * x * x * exp(-0.5*x);
  OCE1PsiOne(oce, c, 0, 0, x, &y_calc);
  ASSERT_NEAR(-0.125,
	      PetscRealPart(k), 0.0001);
  ASSERT_NEAR(PetscRealPart(y_ref),
	      PetscRealPart(y_calc), 0.0001);  
  ASSERT_NEAR(0.0, PetscImaginaryPart(y_calc), 0.0000001);

  
  
  OCE1Destroy(&oce); PFDestroy(&pot);
  MatDestroy(&H); MatDestroy(&V);
  VecDestroy(&c);
  EEPSDestroy(&eps);
}
TEST(TestOCE1DVR, H2plus) {

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
  KSPDestroy(&ksp);  PFDestroy(&slater); 
  VecDestroy(&c_guess); 
  EEPSDestroy(&eps);  

}

TEST(TestOCE1DVR, H_DipOld) {
  /*
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;

  // -- numerical basis --
  Y1s y1s_0, y1s_1;
  ierr = Y1sCreate(comm, &y1s_0); ASSERT_EQ(0, ierr);  
  ierr = Y1sCreate(comm, &y1s_1); ASSERT_EQ(0, ierr);
  ierr = Y1sSet(y1s_0, 0, GERADE, 0);   ASSERT_EQ(0, ierr);
  ierr = Y1sSet(y1s_1, 0, UNGERADE, 1); ASSERT_EQ(0, ierr);

  BPS bps;
  ierr = BPSCreate(comm, &bps); ASSERT_EQ(0, ierr);
  ierr = BPSSetLine(bps, 20.0, 21); ASSERT_EQ(0, ierr);
  DVR dvr;
  ierr = DVRCreate(comm, &dvr); ASSERT_EQ(0, ierr);
  ierr = DVRSetKnots(dvr, 5, bps); ASSERT_EQ(0, ierr);
  ierr = DVRSetUp(dvr); ASSERT_EQ(0, ierr);
  FEMInf fem;
  ierr = FEMInfCreate(comm, &fem); ASSERT_EQ(0, ierr);
  ierr = FEMInfSetDVR(fem, dvr); ASSERT_EQ(0, ierr);
  OCE1 oce_0, oce_1;
  ierr = OCE1Create(comm, &oce_0); ASSERT_EQ(0, ierr);
  ierr = OCE1Create(comm, &oce_1); ASSERT_EQ(0, ierr);
  ierr = OCE1Set(oce_0, fem, y1s_0); ASSERT_EQ(0, ierr);
  ierr = OCE1Set(oce_1, fem, y1s_1); ASSERT_EQ(0, ierr);

  // -- initial state --
  Pot pot; PotCreate(comm, &pot); 
  PotSetCoulombNE(pot, 0, 0.0, -1.0);
  Mat H0; OCE1TMat(oce_0, MAT_INITIAL_MATRIX, &H0);
  Mat V0; OCE1PotMat(oce_0, ROT_SCALAR, pot, MAT_INITIAL_MATRIX, &V0);  
  MatAXPY(H0, 1.0, V0, DIFFERENT_NONZERO_PATTERN);
  EEPS eps0;
  ierr = EEPSCreate(comm, &eps0); ASSERT_EQ(0, ierr);
  ierr = EEPSSetOperators(eps0, H0, NULL); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps0, -0.6); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps0); ASSERT_EQ(0, ierr);
  PetscScalar ene0;
  Vec c0;
  MatCreateVecs(H0, &c0, NULL);
  EEPSGetEigenpair(eps0, 0, &ene0, c0);
  PetscScalar x = 2.5;
  PetscScalar y_ref = 2.0 * x * exp(-x);
  PetscScalar y_calc;
  OCE1PsiOne(oce_0, c0, 0, 0, x, &y_calc);
  ASSERT_NEAR(PetscRealPart(y_ref),
	      PetscRealPart(y_calc), 0.0001);  

  // -- final state --
  Mat H1; OCE1TMat(oce_1, MAT_INITIAL_MATRIX, &H1);
  Mat V1; OCE1PotMat(oce_1, ROT_SCALAR, pot, MAT_INITIAL_MATRIX, &V1);
  MatAXPY(H1, 1.0, V1, DIFFERENT_NONZERO_PATTERN);
  EEPS eps1;
  ierr = EEPSCreate(comm, &eps1); ASSERT_EQ(0, ierr);
  ierr = EEPSSetOperators(eps1, H1, NULL); ASSERT_EQ(0, ierr);
  ierr = EEPSSetTarget(eps1, -0.2); ASSERT_EQ(0, ierr);
  ierr = EEPSSolve(eps1); ASSERT_EQ(0, ierr);
  PetscScalar ene1;
  Vec c1;
  ierr = MatCreateVecs(H1, &c1, NULL); ASSERT_EQ(0, ierr);
  ierr = EEPSGetEigenpair(eps1, 0, &ene1, c1); ASSERT_EQ(0, ierr);
  y_ref = 1.0/sqrt(24.0) * x * x * exp(-0.5*x);
  OCE1PsiOne(oce_1, c1, 0, 0, x, &y_calc);
  ASSERT_NEAR(-0.125,
	      PetscRealPart(ene1), 0.0001);        
  ASSERT_NEAR(PetscRealPart(y_ref),
	      PetscRealPart(y_calc), 0.0001);

  // -- dipole (length) --
  Mat Z;
  ierr = OCE1ZMat(oce_1, oce_0, MAT_INITIAL_MATRIX, &Z); ASSERT_EQ(0, ierr);
  Vec Zc0;
  ierr = MatCreateVecs(Z, NULL, &Zc0); ASSERT_EQ(0, ierr);  
  ierr = MatMult(Z, c0, Zc0); ASSERT_EQ(0, ierr);
  PetscScalar zdip;
  ierr = VecTDot(c1, Zc0, &zdip); ASSERT_EQ(0, ierr);
  ierr = VecMatVecMult(c1, Z, c0, &zdip); ASSERT_EQ(0, ierr);

  // -- dipole (velocity) --
  Mat DZ;
  ierr = OCE1DZMat(oce_1, oce_0, MAT_INITIAL_MATRIX, &DZ); ASSERT_EQ(0, ierr);
  Vec DZc0;
  ierr = MatCreateVecs(DZ, NULL, &DZc0); ASSERT_EQ(0, ierr);
  ierr = MatMult(DZ, c0, DZc0); ASSERT_EQ(0, ierr);
  PetscScalar zdip_v;
  ierr = VecTDot(c1, DZc0, &zdip_v); ASSERT_EQ(0, ierr);
  ierr = VecMatVecMult(c1, DZ, c0, &zdip); ASSERT_EQ(0, ierr);

  // -- see support/hatom_dip.py --
  PetscReal ref = 128.0*sqrt(6)/243.0 * Y1ElePq(1,1,0,0,0);
  ASSERT_NEAR(ref, PetscRealPart(zdip), 0.002);
  ASSERT_NEAR(0.0, PetscImaginaryPart(zdip), 0.002);

  ref = -8.0 * sqrt(6.0) / 81.0 * Y1ElePq(1,1,0,
					  0,  0) ;
  ASSERT_NEAR(ref, PetscRealPart(zdip_v), 0.001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(zdip_v), 0.001);
  
  // -- Finalize --
  FEMInfDestroy(&fem);
  oce_1->fem = NULL; OCE1Destroy(&oce_1);
  oce_0->fem = NULL; OCE1Destroy(&oce_0);
  PFDestroy(&pot);  
  MatDestroy(&H0); MatDestroy(&V0);  EEPSDestroy(&eps0); VecDestroy(&c0);
  MatDestroy(&H1); MatDestroy(&V1);  EEPSDestroy(&eps1); VecDestroy(&c1);
  MatDestroy(&Z);  VecDestroy(&Zc0); MatDestroy(&DZ);    VecDestroy(&DZc0);
  */
  
}
TEST(TestOCE1DVR, H_Dip) {

  // 2p->ks channel
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  
  // -- numerical basis --
  Y1s y1s_0; Y1sCreate(comm, &y1s_0);
  ierr=Y1sSetOne(y1s_0, 0, 0);ASSERT_EQ(0,ierr);
  Y1s y1s_1; Y1sCreate(comm, &y1s_1);
  ierr=Y1sSetOne(y1s_1, 1, 0);ASSERT_EQ(0,ierr);
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps); DVRSetUp(dvr);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);
  OCE1 oce_1s; OCE1Create(comm, &oce_1s); OCE1Set(oce_1s, fem, y1s_0);
  OCE1 oce_2p; OCE1Create(comm, &oce_2p); OCE1Set(oce_2p, fem, y1s_1);

  // -- Fit target --
  // 1s = 2 r exp(-r)
  Pot h_1s; PotCreate(comm, &h_1s);
  PotSetSlater(h_1s, 2.0, 1, 1.0);
  Vec c_1s; OCE1CreateVec(oce_1s, &c_1s);
  ierr = OCE1Fit(oce_1s, h_1s, 0, NULL, c_1s); ASSERT_EQ(0, ierr);

  // 2p = 1/(2sqrt(6)) r**2 exp(-r/2)
  PetscScalar a_2p = 1.0/(2.0 * sqrt(6.0));
  Pot h_2p; PotCreate(comm, &h_2p);
  PotSetSlater(h_2p, a_2p, 2, 0.5);
  Vec c_2p; OCE1CreateVec(oce_2p, &c_2p);
  ierr = OCE1Fit(oce_2p, h_2p, 1, NULL, c_2p);  ASSERT_EQ(0, ierr);

  // compute matrix element
  Mat Z;
  ierr = OCE1ZMat(oce_1s, oce_2p, MAT_INITIAL_MATRIX, &Z); ASSERT_EQ(0, ierr);
  PetscScalar zdip;
  ierr = VecMatVecMult(c_1s, Z, c_2p, &zdip); ASSERT_EQ(0, ierr);

  // -- dipole (velocity) --
  Mat DZ_1s_2p;
  ierr = OCE1DZMat(oce_1s, oce_2p, MAT_INITIAL_MATRIX, &DZ_1s_2p); ASSERT_EQ(0, ierr);
  Mat DZ_2p_1s;
  ierr = OCE1DZMat(oce_2p, oce_1s, MAT_INITIAL_MATRIX, &DZ_2p_1s); ASSERT_EQ(0, ierr);

  PetscScalar zdip_1s_2p;
  VecMatVecMult(c_1s, DZ_1s_2p, c_2p, &zdip_1s_2p);
  PetscScalar zdip_2p_1s;
  VecMatVecMult(c_2p, DZ_2p_1s, c_1s, &zdip_2p_1s);
  // -- see support/hatom_dip.py --
  PetscReal ref_l = 128*sqrt(6)/243 * Y1ElePq(1,1,0,0,0);
  ASSERT_NEAR(ref_l, PetscRealPart(zdip), 0.002);
  ASSERT_NEAR(0.0,   PetscImaginaryPart(zdip), 0.002);

  PetscReal ref_v = 16.0 * sqrt(6.0) / 81.0 * Y1ElePq(1,1,0,
						      0,  0) ;
  ASSERT_NEAR(-ref_v, PetscRealPart(zdip_2p_1s), 0.001);
  ASSERT_NEAR(ref_v, PetscRealPart(zdip_1s_2p), 0.001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(zdip_1s_2p), 0.001);
  ASSERT_NEAR(0.0, PetscImaginaryPart(zdip_2p_1s), 0.001);

  // -- Finalize --
  FEMInfDestroy(&fem);
  oce_1s->fem = NULL; OCE1Destroy(&oce_1s);
  oce_2p->fem = NULL; OCE1Destroy(&oce_2p);
  PFDestroy(&h_1s); VecDestroy(&c_1s);
  PFDestroy(&h_2p); VecDestroy(&c_2p);
  MatDestroy(&Z);  MatDestroy(&DZ_1s_2p); MatDestroy(&DZ_2p_1s);
  
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
