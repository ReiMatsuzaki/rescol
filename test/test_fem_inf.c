#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/mat.h>
#include <rescol/fem_inf.h>
#include <rescol/viewerfunc.h>
#include <rescol/eeps.h>

static char help[] = "Unit test for fem_inf.c \n";

int test1() {

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);

  if(getenv("SHOW_DEBUG")) {
    PetscErrorCode ierr;
    ierr = FEMInfView(fem, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  }

  FEMInfDestroy(&fem);

  return 0;
}
int testH_BSS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 101);
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

  Mat H; FEMInfCreateMat(fem, 1, &H); FEMInfD2R1Mat(fem, H); MatScale(H, -0.5);
  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V); 
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

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
  Vec cs;
  EPSGetConverged(eps->eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  MatCreateVecs(H, &cs, NULL);
  EPSGetEigenpair(eps->eps, 0, &kr, NULL, cs, NULL);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));
  
  FEMInfDestroy(&fem);
  MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  EEPSDestroy(&eps); VecDestroy(&x0[0]); VecDestroy(&cs);
  return 0;
}
int testH_BSS_accurate() {
  // see    HBauchau, E.Comrmier, JDecleva et al.
  //           Reports on Progress in Physics 64, 1815
  // see 2016/4/20

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

  return 0;
}
int testH_PI_BSS() {

  // from calc/stoh/1skp/l_5
  //1.88562800720386      , -0.362705406693342

  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_PI", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 101);
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
  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V); 
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
  ierr = VecSetType(m, "seq"); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c); CHKERRQ(ierr);

  Vec s_m; VecDuplicate(m, &s_m); 
  ierr = MatMult(S, m, s_m); CHKERRQ(ierr);  

  PetscScalar alpha;
  ierr = VecTDot(c, s_m, &alpha); CHKERRQ(ierr);
  ASSERT_DOUBLE_EQ(1.88562800720386, creal(alpha));
  ASSERT_DOUBLE_EQ(-0.362705406693342, cimag(alpha));

  return 0;

}
int testPOT_BSS() {
  
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "POT_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 20.0, 21, 5.0);
  BSS bss; BSSCreate(comm, &bss); 
  BSSSetKnots(bss, 5, bps); BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

  Pot r2inv; PotCreate(comm, &r2inv); PotSetPower(r2inv, 1.0, -2.0); 
  Mat A; FEMInfCreateMat(fem, 1, &A); FEMInfPotR1Mat(fem, r2inv, A);
  Mat B; FEMInfCreateMat(fem, 1, &B); FEMInfR2invR1Mat(fem, B);

  MatAXPY(A, -1.0, B, DIFFERENT_NONZERO_PATTERN);

  PetscReal norm;
  MatNorm(A, NORM_1, &norm);
  ASSERT_DOUBLE_NEAR(norm, 0.0, 0.00000001);

  FEMInfDestroy(&fem); PFDestroy(&r2inv); MatDestroy(&A); MatDestroy(&B);

  return 0;
}
int testH_DVR() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_DVR", NULL);
  
  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 20.0, 21, 5.0);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, 5, bps);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetDVR(fem, dvr);

  Mat H; FEMInfCreateMat(fem, 1, &H); FEMInfD2R1Mat(fem, H); MatScale(H, -0.5);
  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V); 
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  EPS eps; 
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, H, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSJD);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  
  EPSSetTarget(eps, -0.6);
  EPSSolve(eps);
  
  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  FEMInfDestroy(&fem);
  MatDestroy(&H); MatDestroy(&V);
  EPSDestroy(&eps);
  
  return 0;
}
int testFit_BSS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "FIT_BSS", NULL);

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 301);
  CScaling c_scaling; CScalingCreate(comm, &c_scaling);
  CScalingSetSharpECS(c_scaling, 70.0, 20.0);
  
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, 5, bps);
  BSSSetUp(bss);
  BSSSetCScaling(bss, c_scaling); BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);

  Pot sto; PotCreate(comm, &sto); PotSetSlater(sto, 1.1, 2, 1.2);

  KSP ksp; KSPCreate(comm, &ksp); 

  Vec c;   VecCreate(comm, &c); VecSetType(c, "seq");

  FEMInfFit(fem, sto, ksp, c);
  PetscScalar x = 2.2;
  PetscScalar y_calc;
  FEMInfPsi(fem, c, x, &y_calc);

  PetscScalar xs[1] = {x};
  PetscScalar ys[1];
  PFApply(sto, 1, xs, ys);
  PetscScalar y_ref = ys[0];
  
  ASSERT_DOUBLE_NEAR(y_ref, y_calc, pow(10.0, -6.0));

  PFDestroy(&sto);
  VecDestroy(&c);
  KSPDestroy(&ksp);
  FEMInfDestroy(&fem);
  
  
  return 0;
}


int main(int argc, char **args) {
  
  
  SlepcInitialize(&argc, &args, (char*)0, help);
  PetscErrorCode ierr;

  test1();    
  testH_BSS();
  // PetscErrorCode ierr;
  //  ierr = testH_BSS_accurate(); CHKERRQ(ierr);
  testPOT_BSS();
  testH_DVR();
  testFit_BSS();
  ierr = testH_PI_BSS(); CHKERRQ(ierr);

  SlepcFinalize();
  return 0;

}

