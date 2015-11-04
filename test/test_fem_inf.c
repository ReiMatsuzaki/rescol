#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/mat.h>
#include <rescol/fem_inf.h>

static char help[] = "Unit test for fem_inf.c \n";

int test1() {

  PetscErrorCode ierr;
  
  printf("AAA\n");
  
  BPS bps;
  BPSCreate(&bps, PETSC_COMM_SELF); BPSSetLine(bps, 5.0, 6);
  DVR dvr;
  DVRCreate(&dvr, 4, bps, PETSC_COMM_SELF);
  FEMInf fem;
  FEMInfCreateDVR(&fem, dvr);

  ierr = FEMInfFPrintf(fem, stdout, 0); CHKERRQ(ierr);


  return 0;
}
int testH_BSS() {

  BPS bps; BPSCreate(&bps, PETSC_COMM_SELF); BPSSetExp(bps, 20.0, 21, 5.0);
  BSS bss; BSSCreate(&bss, 5, bps, NULL, PETSC_COMM_SELF);
  FEMInf fem; FEMInfCreateBSS(&fem, bss);

  printf("\n");
  FEMInfFPrintf(fem, stdout, 0);
  printf("\n");

  Mat H;
  FEMInfSetD2R1Mat(fem, &H); MatScale(H, -0.5);

  Mat V;
  FEMInfSetENR1Mat(fem, 0, 0.0, &V);
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  Mat S;
  FEMInfSetSR1Mat(fem, &S); 

  EPS eps; 
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, H, S);
  EPSSetProblemType(eps, EPS_GHEP);
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

  

  return 0;
}
int testH_DVR() {

  BPS bps;
  BPSCreate(&bps, PETSC_COMM_SELF); BPSSetExp(bps, 20.0, 21, 5.0);
  DVR dvr;
  DVRCreate(&dvr, 5, bps, PETSC_COMM_SELF);
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
  testH_DVR();
  SlepcFinalize();
  return 0;
}

