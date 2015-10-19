#include "mat.h"
#include "bspline.h"

static char help[] = "solve hydrogen atom eigenvalue problem";

/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_mpd : maximum dim of projected problem
  
  -bss_order : 
  -bss_rmax  : 
  -bss_knots_num :  
*/

int main(int argc, char **args) {

  PetscErrorCode ierr;
  BSS bss;
  MPI_Comm comm = PETSC_COMM_SELF;
  char target_dir[100]  = ".";
  PetscInt L = 0;
  
  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(comm, "", "solve h atom eigen problem", "none");
  PrintTimeStamp(comm, "Init1", NULL);
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PrintTimeStamp(comm, "Init2", NULL);
  PetscOptionsGetInt(NULL, "-L", &L, NULL);
  PrintTimeStamp(comm, "Init3", NULL);
  ierr = BSSCreateFromOptions(&bss, comm);  CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init4", NULL);
  PetscOptionsEnd();

  // Matrix
  Mat H, S, tmp;
  PrintTimeStamp(comm, "Mat", NULL);
  BSSInitR1Mat(bss, PETSC_COMM_SELF, &H);
  BSSInitR1Mat(bss, PETSC_COMM_SELF, &S);

  BSSInitR1Mat(bss, PETSC_COMM_SELF, &H);
  BSSCalcD2R1Mat(bss, H, INSERT_VALUES);
  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  MatScale(H, -0.5);
    
  BSSInitR1Mat(bss, PETSC_COMM_SELF, &tmp);
  BSSCalcENR1Mat(bss, 0, 0.0, tmp, INSERT_VALUES);
  MatAssemblyBegin(tmp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tmp, MAT_FINAL_ASSEMBLY);
  MatAXPY(H, -1.0, tmp, SAME_NONZERO_PATTERN);
  MatDestroy(&tmp);

  BSSCalcSR1Mat(bss, S, INSERT_VALUES);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);  

  // Solve
  EPS eps; 
  PrintTimeStamp(comm, "EPS", NULL);
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, H, S);
  EPSSetProblemType(eps, EPS_GHEP);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6);
  EPSSolve(eps);  

  // Output
  PetscPrintf(PETSC_COMM_SELF, "START H atom\n");
  BSSFPrintf(bss, PETSC_COMM_SELF, stdout, 0);

  int nconv;
  PetscScalar kr, ki;
  Vec xr, xi;
  MatCreateVecs(H, NULL, &xr); MatCreateVecs(H, NULL, &xi);
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    PetscPrintf(PETSC_COMM_SELF, "eig%i: %f\n", i, kr);
  }

  // Finalize
  VecDestroy(&xi); VecDestroy(&xr);
  MatDestroy(&H); MatDestroy(&S);
  EPSDestroy(&eps);
  BSSDestroy(&bss);
  SlepcFinalize();
  return 0;
}
