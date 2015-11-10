#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/mat.h>
#include <rescol/fem_inf.h>

static char help[] = "Unit test for fem_inf.c \n";

int test1() {

  BPS bps;
  BPSCreate(&bps, PETSC_COMM_SELF); BPSSetLine(bps, 5.0, 6);
  DVR dvr;
  DVRCreate(&dvr, 4, bps, PETSC_COMM_SELF);
  FEMInf fem;
  FEMInfCreateDVR(&fem, dvr);

  if(getenv("SHOW_DEBUG")) {
    PetscErrorCode ierr;
    ierr = FEMInfFPrintf(fem, stdout, 0); CHKERRQ(ierr);
  }

  return 0;
}
int testH_BSS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "H_BSS", NULL);

  BPS bps; BPSCreate(&bps, PETSC_COMM_SELF); BPSSetLine(bps, 100.0, 101);
  Scaler scaler; ScalerCreateSharpECS(&scaler, comm, 70.0, 20.0);
  BSS bss; BSSCreate(&bss, 5, bps, scaler, comm);
  FEMInf fem; FEMInfCreateBSS(&fem, bss);

  if(getenv("SHOW_DEBUG")) {
    printf("\n");
    printf("SHOW_DEBUG = %s\n", getenv("SHOW_DEBUG"));
    FEMInfFPrintf(fem, stdout, 0);
    printf("\n");
  }

  Mat H;
  FEMInfSetD2R1Mat(fem, &H); MatScale(H, -0.5);

  Mat V;
  FEMInfSetENR1Mat(fem, 0, 0.0, &V);
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  Mat S;
  FEMInfSetSR1Mat(fem, &S); 

  EPS eps; 
  EPSCreateForBoundState(&eps, comm, H, S, -0.6);
  Vec x0[1]; MatCreateVecs(H, &x0[0], NULL); 
  int num; FEMInfGetSize(fem, &num);
  for(int i = 0; i < num; i++) {
    if(i < 400)
      VecSetValue(x0[0], i, 0.5, INSERT_VALUES);
  }
  VecAssemblyBegin(x0[0]); VecAssemblyEnd(x0[0]);
  EPSSetInitialSpace(eps, 1, x0);
  EPSSolve(eps);
  
  int nconv;
  PetscScalar kr;
  Vec cs;
  EPSGetConverged(eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  MatCreateVecs(H, &cs, NULL);
  EPSGetEigenpair(eps, 0, &kr, NULL, cs, NULL);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  return 0;
}
int testPOT_BSS() {
  
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "POT_BSS", NULL);
  BPS bps; BPSCreate(&bps, comm); BPSSetExp(bps, 20.0, 21, 5.0);
  BSS bss; BSSCreate(&bss, 5, bps, NULL, comm);
  FEMInf fem; FEMInfCreateBSS(&fem, bss);

  POT r2inv; POTPowerCreate(&r2inv, 1.0, -2.0);
  Mat A; FEMInfSetPOTR1Mat(fem, r2inv, &A);
  Mat B; FEMInfSetR2invR1Mat(fem, &B);

  MatAXPY(A, -1.0, B, DIFFERENT_NONZERO_PATTERN);

  PetscReal norm;
  MatNorm(A, NORM_1, &norm);
  ASSERT_DOUBLE_NEAR(norm, 0.0, 0.00000001);
  
  return 0;
}
int testH_DVR() {

  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "POT_BSS", NULL);
  
  BPS bps;
  BPSCreate(&bps, comm); BPSSetExp(bps, 20.0, 21, 5.0);
  DVR dvr;
  DVRCreate(&dvr, 5, bps, comm);
  FEMInf fem;
  FEMInfCreateDVR(&fem, dvr);

  Mat H;
  FEMInfSetD2R1Mat(fem, &H); MatScale(H, -0.5);

  Mat V;
  FEMInfSetENR1Mat(fem, 0, 0.0, &V);
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

