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

int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);

  test1();    
  testH_BSS();
  testPOT_BSS();
  testH_DVR();

  SlepcFinalize();
  return 0;

}

