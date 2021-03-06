#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include "../include/math.h"
#include "../include/mat.h"
#include "../include/fem_inf.h"
#include "../include/viewerfunc.h"
#include "../include/eeps.h"

static char help[] = "Unit test for fem_inf.c \n";

int testH_BSS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 201);
  CScaling c_scaling; 
  CScalingCreate(comm, &c_scaling); CScalingSetSharpECS(c_scaling, 70.0, 20.0);
  BSS bss; BSSCreate(comm, &bss); 
  BSSSetKnots(bss, 5, bps); BSSSetCScaling(bss, c_scaling); BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

  if(getenv("SHOW_DEBUG")) {
    printf("\n");
    printf("SHOW_DEBUG = %s\n", getenv("SHOW_DEBUG"));
    FEMInfView(fem, PETSC_VIEWER_STDOUT_SELF);
    printf("\n");
  }

  Pot pot; PotCreate(comm, &pot);
  PotSetCoulombNE(pot, 0, 0.0, -1.0);

  Mat H; FEMInfCreateMat(fem, 1, &H); FEMInfD2R1Mat(fem, H); MatScale(H, -0.5);
  //  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V);
  Mat V; FEMInfCreateMat(fem, 1, &V);
  FEMInfPotR1Mat(fem, pot, V);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);

  Mat S; FEMInfCreateMat(fem, 1, &S); FEMInfSR1Mat(fem, S);

  EEPS eps; EEPSCreate(comm, &eps);
  EEPSSetOperators(eps, H, S);
  EEPSSetTarget(eps, -0.6);

  Vec x0[1]; MatCreateVecs(H, &x0[0], NULL); 
  int num; FEMInfGetSize(fem, &num);
  for(int i = 0; i < num; i++) {
    if(i < 400)
      VecSetValue(x0[0], i, 0.5, INSERT_VALUES);
  }
  VecAssemblyBegin(x0[0]); VecAssemblyEnd(x0[0]);
  EPSSetInitialSpace(eps->eps, 1, x0);
  EEPSSolve(eps);
  
  int nconv;
  PetscScalar kr;
  
  EPSGetConverged(eps->eps, &nconv);
  ASSERT_TRUE(nconv > 0);  
  EPSGetEigenpair(eps->eps, 0, &kr, NULL, NULL, NULL);
  
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  Vec cs;
  MatCreateVecs(H, &cs, NULL);
  EEPSGetEigenvector(eps, 0, cs);

  PetscReal x = 1.1;
  PetscScalar xc = 1.1;
  PetscScalar y, dy;
  FEMInfPsiOne(fem, cs, xc, &y); FEMInfDerivPsiOne(fem, cs, x, &dy);
  ASSERT_DOUBLE_NEAR(2.0*x*exp(-x), creal(y), pow(10.0, -5.0));
  ASSERT_DOUBLE_NEAR(2.0*(1.0-x)*exp(-x), creal(dy), pow(10.0, -3.0));
  
  FEMInfDestroy(&fem);
  PFDestroy(&pot);
  MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  EEPSDestroy(&eps); VecDestroy(&x0[0]); VecDestroy(&cs);
  return 0;
}
int testH_BSS_accurate() {
  // see    HBauchau, E.Comrmier, JDecleva et al.
  //           Reports on Progress in Physics 64, 1815
  // see 2016/4/20

  /*
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_BSS2", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 100.0, 101, 5.0);
  CScaling c_scaling; 
  ierr = CScalingCreate(comm, &c_scaling);           CHKERRQ(ierr);
  ierr = CScalingSetSharpECS(c_scaling, 70.0, 20.0); CHKERRQ(ierr);
  BSS bss; BSSCreate(comm, &bss); 
  ierr = BSSSetKnots(bss, 9, bps); CHKERRQ(ierr);
  BSSSetCScaling(bss, c_scaling); BSSSetUp(bss); CHKERRQ(ierr);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);
  
  // compute matrix
  printf("a\n");
  Mat H;
  ierr = FEMInfCreateMat(fem, 1, &H);  CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(fem, H); CHKERRQ(ierr);
  ierr = MatScale(H, -0.5); CHKERRQ(ierr);
  Mat V;
  ierr = FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V); CHKERRQ(ierr);
  ierr = MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  Mat S;
  ierr = FEMInfCreateMat(fem, 1, &S); CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(fem, S); CHKERRQ(ierr);
  
  // create initial space
  printf("b\n");
  
  Pot psi0; PotCreate(comm, &psi0); PotSetSlater(psi0, 2.0, 1, 1.0);
  int n_init_space = 1;

  printf("c\n");
  Vec *xs;
  KSP ksp; KSPCreate(comm, &ksp); KSPSetFromOptions(ksp);
  ierr = PetscMalloc1(n_init_space, &xs); CHKERRQ(ierr);
  ierr = VecCreate(comm, &xs[0]); VecSetType(xs[0], "seq");
  ierr = FEMInfFit(fem, psi0, ksp, xs[0]); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  printf("d\n");
  EEPS eps; EEPSCreate(comm, &eps);
  ierr = EEPSSetOperators(eps, H, S); CHKERRQ(ierr);
  ierr = EEPSSetTarget(eps, -0.6); CHKERRQ(ierr);
  ierr = EPSSetType(eps->eps, EPSJD); CHKERRQ(ierr);
  //ierr = EPSSetType(eps->eps, EPSARNOLDI); CHKERRQ(ierr);
  //ierr = EPSSetType(eps->eps, EPSGD); CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(eps->eps, n_init_space, xs); CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps->eps, 4, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  //ierr = EPSSetTolerances(eps->eps, pow(10.0, -10.0), 100); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps->eps); CHKERRQ(ierr);

  printf("f\n");
  ierr = EPSSolve(eps->eps); CHKERRQ(ierr);

  int nconv;
  PetscScalar kr;
  Vec cs;
  ierr = EPSGetConverged(eps->eps, &nconv); CHKERRQ(ierr);
  printf("nconv = %d\n", nconv);
  ASSERT_TRUE(nconv > 5);
  
  
  EPSGetEigenpair(eps->eps, 0, &kr, NULL, NULL, NULL); 
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -9.0));
  
  EPSGetEigenpair(eps->eps, 1, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.125, kr, pow(10.0, -9.0));

  EPSGetEigenpair(eps->eps, 2, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-1.0/(2.0*3.0*3.0), kr, pow(10.0, -9.0));

  EPSGetEigenpair(eps->eps, 3, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-1.0/(2.0*4.0*4.0), kr, pow(10.0, -9.0));

  EPSGetEigenpair(eps->eps, 4, &kr, NULL, NULL, NULL);
  //ASSERT_DOUBLE_NEAR(-1.0(2.0*5.0*5.0), kr, pow(10.0, -7.0));
  ASSERT_DOUBLE_NEAR(-0.0199999715, kr, pow(10.0, -10.0));

  
  FEMInfDestroy(&fem);
  MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  EEPSDestroy(&eps);  VecDestroy(&cs);
  // VecDestroy(&xs[0]);
  */
  
  return 0;
}
int testH_PI_BSS() {

  // from calc/stoh/1skp/l_5
  //1.88562800720386      , -0.362705406693342

  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_PI_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 201);
  CScaling c_scaling; 
  CScalingCreate(comm, &c_scaling); CScalingSetSharpECS(c_scaling, 70.0, 20.0);
  BSS bss; BSSCreate(comm, &bss); 
  BSSSetKnots(bss, 5, bps); BSSSetCScaling(bss, c_scaling); BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

  if(getenv("SHOW_DEBUG")) {
    printf("\n");
    printf("SHOW_DEBUG = %s\n", getenv("SHOW_DEBUG"));
    FEMInfView(fem, PETSC_VIEWER_STDOUT_SELF);
    printf("\n");
  }

  Mat L; FEMInfCreateMat(fem, 1, &L); FEMInfD2R1Mat(fem, L);
  MatScale(L, -0.5);
  Pot pot; PotCreate(comm, &pot);
  PotSetCoulombNE(pot, 0, 0.0, -1.0);  
  Mat V; FEMInfCreateMat(fem, 1, &V); //FEMInfENR1Mat(fem, 0, 0.0, V);
  FEMInfPotR1Mat(fem, pot, V);
  MatAXPY(L, -1.0, V, DIFFERENT_NONZERO_PATTERN);
  Pot r2; PotCreate(comm, &r2); PotSetPower(r2, 1.0, -2);
  Mat LV; FEMInfCreateMat(fem, 1, &LV); FEMInfPotR1Mat(fem, r2, LV);
  //  Mat LV; FEMInfCreateMat(fem, 1, &LV); FEMInfR2invR1Mat(fem, LV); 
  MatAXPY(L, 1.0, LV, DIFFERENT_NONZERO_PATTERN);
  Mat S; FEMInfCreateMat(fem, 1, &S); FEMInfSR1Mat(fem, S);
  MatAXPY(L, -0.5, S, DIFFERENT_NONZERO_PATTERN);

  Pot driv; PotCreate(comm, &driv); PotSetSlater(driv, 2.0, 2, 1.0);
  Vec m; FEMInfCreateVec(fem, 1, &m); FEMInfPotR1Vec(fem, driv, m);
  
  KSP ksp;
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, L, L); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  Vec c;
  ierr = VecCreate(comm, &c); CHKERRQ(ierr);
  int n; FEMInfGetSize(fem, &n); 
  ierr = VecSetSizes(c, n, n); CHKERRQ(ierr);
  ierr = VecSetType(c, "seq"); CHKERRQ(ierr);
  //ierr = VecSetType(m, "seq"); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c); CHKERRQ(ierr);

  PetscScalar alpha;
  ierr = VecTDot(c, m, &alpha); CHKERRQ(ierr);
  
  PetscScalar a_ref = -5.65688402161+I*1.08811622008;
  ASSERT_SCALAR_EQ(a_ref, alpha);
  return 0;

}
int testH_PI_DVR() {

  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_PI_DVR", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 201);
  CScaling cscaling; CScalingCreate(comm, &cscaling); 
  CScalingSetSharpECS(cscaling, 70.0, 20.0*M_PI/180.0);
  DVR dvr; DVRCreate(comm, &dvr);
  DVRSetKnots(dvr, 5, bps); 
  DVRSetCScaling(dvr, cscaling);
  DVRSetUp(dvr);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);

  if(getenv("SHOW_DEBUG")) {
    printf("\n");
    printf("SHOW_DEBUG = %s\n", getenv("SHOW_DEBUG"));
    FEMInfView(fem, PETSC_VIEWER_STDOUT_SELF);
    printf("\n");
  }

  Mat L; FEMInfCreateMat(fem, 1, &L); FEMInfD2R1Mat(fem, L);
  MatScale(L, -0.5);
  
  //Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V);
  Pot pot; PotCreate(comm, &pot);
  PotSetCoulombNE(pot, 0, 0.0, -1.0);
  Mat V; FEMInfCreateMat(fem, 1, &V);
  FEMInfPotR1Mat(fem, pot, V);
  MatAXPY(L, 1.0, V, DIFFERENT_NONZERO_PATTERN);
  
  Pot r2; PotCreate(comm, &r2); PotSetPower(r2, 1.0, -2);
  Mat LV; FEMInfCreateMat(fem, 1, &LV); FEMInfPotR1Mat(fem, r2, LV);
  //  Mat LV; FEMInfCreateMat(fem, 1, &LV); FEMInfR2invR1Mat(fem, LV); 
  MatAXPY(L, 1.0, LV, DIFFERENT_NONZERO_PATTERN);
  Mat S; FEMInfCreateMat(fem, 1, &S); FEMInfSR1Mat(fem, S);
  MatAXPY(L, -0.5, S, DIFFERENT_NONZERO_PATTERN);

  Pot driv; PotCreate(comm, &driv); PotSetSlater(driv, 2.0, 2, 1.0);
  Vec m; FEMInfCreateVec(fem, 1, &m); FEMInfPotR1Vec(fem, driv, m);
  
  KSP ksp;
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, L, L); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  Vec c;
  ierr = VecCreate(comm, &c); CHKERRQ(ierr);
  int n; FEMInfGetSize(fem, &n); 
  ierr = VecSetSizes(c, n, n); CHKERRQ(ierr);
  ierr = VecSetType(c, "seq"); CHKERRQ(ierr);
  //ierr = VecSetType(m, "seq"); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c); CHKERRQ(ierr);

  PetscScalar alpha;
  ierr = VecTDot(c, m, &alpha); CHKERRQ(ierr);

  ASSERT_DOUBLE_NEAR(-5.65688402161, creal(alpha), 0.000001);
  ASSERT_DOUBLE_NEAR(1.08811622008, cimag(alpha),  0.000001);

  FEMInfDestroy(&fem);
  MatDestroy(&L); PFDestroy(&pot); 
  MatDestroy(&V); MatDestroy(&LV); MatDestroy(&S);
  PFDestroy(&r2); PFDestroy(&driv); VecDestroy(&m); 
  KSPDestroy(&ksp); VecDestroy(&c); 

  return 0;

}
int testH_DVR() {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_DVR", NULL);
  BPS bps;
  ierr = BPSCreate(comm, &bps); CHKERRQ(ierr);
  ierr = BPSSetExp(bps, 20.0, 41, 5.0);CHKERRQ(ierr);
  CScaling cscaling;
  CScalingCreate(comm, &cscaling);
  CScalingSetSharpECS(cscaling, 15.0, 20.0*M_PI/180.0);
  DVR dvr;
  ierr = DVRCreate(comm, &dvr);CHKERRQ(ierr);
  ierr = DVRSetKnots(dvr, 8, bps);CHKERRQ(ierr);
  ierr = DVRSetCScaling(dvr, cscaling);CHKERRQ(ierr);
  ierr = DVRSetUp(dvr); CHKERRQ(ierr);
  FEMInf fem;
  ierr = FEMInfCreate(comm, &fem);CHKERRQ(ierr);
  ierr = FEMInfSetDVR(fem, dvr);CHKERRQ(ierr);

  Mat H;
  ierr = FEMInfCreateMat(fem, 1, &H);CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(fem, H);CHKERRQ(ierr);
  ierr = MatScale(H, -0.5);CHKERRQ(ierr);
  Pot pot; PotCreate(comm, &pot);
  PotSetCoulombNE(pot, 0, 0.0, -1.0);
  Mat V;
  ierr = FEMInfCreateMat(fem, 1, &V);CHKERRQ(ierr);
  //  ierr = FEMInfENR1Mat(fem, 0, 0.0, V); CHKERRQ(ierr);
  //  ierr = MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = FEMInfPotR1Mat(fem, pot, V); CHKERRQ(ierr);
  ierr = MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  KSP ksp;
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  Pot slater;
  ierr = PotCreate(comm, &slater); CHKERRQ(ierr);
  ierr = PotSetSlater(slater, 1.9, 1, 1.0); CHKERRQ(ierr);
  //  Vec c_fit; VecCreate(comm, &c_fit); VecSetType(c_fit, "seq");
  Vec c_fit; FEMInfCreateVec(fem, 1, &c_fit); 
  ierr = FEMInfFit(fem, slater, ksp, c_fit); CHKERRQ(ierr);

  EEPS eps;
  ierr = EEPSCreate(PETSC_COMM_SELF, &eps);CHKERRQ(ierr);
  ierr = EEPSSetOperators(eps, H, NULL); CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(eps->eps, 1, &c_fit); CHKERRQ(ierr);
  ierr = EPSSetType(eps->eps, EPSGD);CHKERRQ(ierr);
  ierr = EEPSSetTarget(eps, -0.6);CHKERRQ(ierr);
  ierr = EEPSSolve(eps);CHKERRQ(ierr);

  int nconv;
  PetscScalar kr;
  ierr = EPSGetConverged(eps->eps, &nconv);CHKERRQ(ierr);
  ASSERT_TRUE(nconv > 0);
  ierr = EPSGetEigenpair(eps->eps, 0, &kr, NULL, NULL, NULL);CHKERRQ(ierr);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  Vec cs;
  ierr = MatCreateVecs(H, &cs, NULL); CHKERRQ(ierr);
  ierr = EEPSGetEigenvector(eps, 0, cs); CHKERRQ(ierr);
  //  VecView(cs, PETSC_VIEWER_STDOUT_SELF);

  PetscReal x = 1.1;
  PetscScalar y, dy;

  ierr = FEMInfPsiOne(fem, cs, x, &y); CHKERRQ(ierr);

  ierr = FEMInfDerivPsiOne(fem, cs, x, &dy); CHKERRQ(ierr);
  ASSERT_SCALAR_NEAR(2.0*x*exp(-x), y, pow(10.0, -6.0)); 
  ASSERT_SCALAR_NEAR(2.0*exp(-x)-2.0*x*exp(-x), dy, pow(10.0, -6.0));

  // -- Destroy --
  FEMInfDestroy(&fem); MatDestroy(&H); PFDestroy(&pot);
  MatDestroy(&V);     
  KSPDestroy(&ksp);    PFDestroy(&slater);
  VecDestroy(&c_fit);  EEPSDestroy(&eps);
  VecDestroy(&cs);    

  return 0;
}
int testFit_BSS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "FIT_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 301);
  //CScaling c_scaling; CScalingCreate(comm, &c_scaling);
  //CScalingSetSharpECS(c_scaling, 70.0, 20.0);
  
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, 5, bps);
  BSSSetUp(bss);
  //  BSSSetCScaling(bss, c_scaling);
  BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

  Pot sto; PotCreate(comm, &sto); PotSetSlater(sto, 1.1, 2, 1.2);

  KSP ksp; KSPCreate(comm, &ksp); 

  Vec c;   FEMInfCreateVec(fem, 1, &c); VecSetType(c, "seq");

  FEMInfFit(fem, sto, ksp, c);
  PetscScalar x = 2.2;
  PetscScalar y_calc;
  FEMInfPsiOne(fem, c, x, &y_calc);

  PetscScalar xs[1] = {x};
  PetscScalar ys[1] = {0.0};
  PFApply(sto, 1, xs, ys);
  PetscScalar y_ref = ys[0];
  
  ASSERT_DOUBLE_NEAR(1.1*x*x*exp(-1.2*x), y_calc, pow(10.0, -6.0));
  ASSERT_DOUBLE_NEAR(y_ref, y_calc, pow(10.0, -6.0));

  FEMInfDestroy(&fem);  
  PFDestroy(&sto);
  KSPDestroy(&ksp);
  VecDestroy(&c);
  
  return 0;
}
int testFit_DVR() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "FIT_DVR", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 301);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps);
  DVRSetUp(dvr);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);

  Pot sto; PotCreate(comm, &sto); PotSetSlater(sto, 1.1, 2, 1.2);
  KSP ksp; KSPCreate(comm, &ksp); 
  Vec c;   FEMInfCreateVec (fem, 1, &c); VecSetType(c, "seq");

  FEMInfFit(fem, sto, ksp, c);
  PetscScalar x = 2.2;
  PetscScalar y_calc;
  FEMInfPsiOne(fem, c, x, &y_calc);

  PetscScalar xs[1] = {x};
  PetscScalar ys[1] = {0.0};
  PFApply(sto, 1, xs, ys);
  PetscScalar y_ref = ys[0];
  
  ASSERT_DOUBLE_NEAR(1.1*x*x*exp(-1.2*x), y_calc, pow(10.0, -6.0));
  ASSERT_DOUBLE_NEAR(y_ref, y_calc, pow(10.0, -6.0));

  // ---- Destroy ----
  FEMInfDestroy(&fem);
  PFDestroy(&sto);
  KSPDestroy(&ksp);
  VecDestroy(&c);
  
  return 0;
}
int testDZMat_DVR() {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "DZMat", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 301);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps); DVRSetUp(dvr);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);

  // see support/hatom_dip.py
  PetscScalar a0 = 2.0;
  PetscScalar a1 = 1.0/(2.0*sqrt(6.0));
  PetscInt    n0 = 1;
  PetscInt    n1 = 2;    
  PetscScalar z0 = 1.0;
  PetscScalar z1 = 0.5;

  Pot sto0; PotCreate(comm, &sto0); PotSetSlater(sto0, a0, n0, z0);
  Pot sto1; PotCreate(comm, &sto1); PotSetSlater(sto1, a1, n1, z1);
  Vec c0; FEMInfCreateVec (fem, 1, &c0); VecSetType(c0, "seq");
  FEMInfFit(fem, sto0, NULL, c0);
  Vec c1; FEMInfCreateVec (fem, 1, &c1); VecSetType(c1, "seq");
  FEMInfFit(fem, sto1, NULL, c1);  

  Mat dz;
  ierr = FEMInfCreateMat(fem, 1, &dz); CHKERRQ(ierr);
  ierr = FEMInfD1R1Mat(fem, dz); CHKERRQ(ierr);
  Vec dz_c0;
  ierr = MatCreateVecs(dz, NULL, &dz_c0); CHKERRQ(ierr);
  ierr = MatMult(dz, c0, dz_c0); CHKERRQ(ierr);
  PetscScalar dip;
  ierr = VecTDot(c1, dz_c0, &dip); CHKERRQ(ierr);
  PetscScalar ref = -8.0 * sqrt(6.0) / 81.0;
  
  ASSERT_SCALAR_NEAR(ref, dip, 0.000001);

  FEMInfDestroy(&fem);
  PFDestroy(&sto0);
  PFDestroy(&sto1);
  VecDestroy(&c0);
  VecDestroy(&c1);
  MatDestroy(&dz);
  VecDestroy(&dz_c0);
  
  return 0;
}
int testCopy_DVR() {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;  
  PrintTimeStamp(comm, "Copy_DVR", NULL);
  
  BPS bps;
  ierr = BPSCreate(comm, &bps); CHKERRQ(ierr);
  ierr = BPSSetExp(bps, 20.0, 41, 5.0);CHKERRQ(ierr);
  CScaling cscaling;
  CScalingCreate(comm, &cscaling);
  CScalingSetSharpECS(cscaling, 15.0, 20.0*M_PI/180.0);
  DVR dvr;
  ierr = DVRCreate(comm, &dvr);CHKERRQ(ierr);
  ierr = DVRSetKnots(dvr, 8, bps);CHKERRQ(ierr);
  ierr = DVRSetCScaling(dvr, cscaling);CHKERRQ(ierr);
  ierr = DVRSetUp(dvr); CHKERRQ(ierr);
  FEMInf fem;
  ierr = FEMInfCreate(comm, &fem);CHKERRQ(ierr);
  ierr = FEMInfSetDVR(fem, dvr);CHKERRQ(ierr);

  FEMInf fem2;
  ierr = FEMInfDuplicate(fem, &fem2); CHKERRQ(ierr);
  ierr = FEMInfCopy(fem, fem2); CHKERRQ(ierr);

  Mat H;
  ierr = FEMInfCreateMat(fem, 1, &H);CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(fem, H);CHKERRQ(ierr);

  Mat H1;
  ierr = FEMInfCreateMat(fem, 1, &H1); CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(fem, H1); CHKERRQ(ierr);

  int i[1] = {0}; 
  int j[1] = {0};
  PetscScalar v[1], v1[1];
  MatGetValues(H, 1, i, 1, j, v);
  MatGetValues(H1,1, i, 1, j, v1);
  ASSERT_SCALAR_EQ(v[0], v1[0]);


  FEMInfDestroy(&fem);
  FEMInfDestroy(&fem2);
  MatDestroy(&H);
  MatDestroy(&H1);
  return 0;

}
int main(int argc, char **args) {
  
  SlepcInitialize(&argc, &args, (char*)0, help);

  
  testH_BSS();
  /// testH_BSS_accurate();
  /// testH_PI_BSS();

  testH_PI_DVR();
  testH_DVR();
  testFit_BSS();
  testFit_DVR();
  testDZMat_DVR();
  testCopy_DVR();
  
  SlepcFinalize();
  return 0;

}

