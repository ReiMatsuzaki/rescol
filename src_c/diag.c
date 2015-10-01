#include "mat.h"

static char help[] = "Read h_mat and s_mat written in COO format and solve eigen value problem";

int main (int argc, char **args)
{

  PetscErrorCode ierr;
  char h_path[256];
  char s_path[256];
  Mat H, S;
  Vec xs, ys;
  EPS eps;

  EPSType type;
  PetscInt its, nev, maxit, nconv, i;
  PetscReal tol, eig ,im_eig, error;

  // Initialization
  ierr = SlepcInitialize(&argc, &args, (char*)0, help);

  // Read from options 
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "test for mat.c options", "none");
  PetscOptionsString("-h-mat", "path for H matrix data file", "", 
		     h_path, h_path, 256, NULL);
  PetscOptionsString("-s-mat", "path for S matrix data file", "", 
		     s_path, s_path, 256, NULL);
  PetscOptionsEnd();

  // Matrix and creation from file
  MatCreateFromCOOFormatFile(h_path, &H);
  MatCreateFromCOOFormatFile(s_path, &S);
  MatView(H, PETSC_VIEWER_STDOUT_SELF);
  ierr = MatCreateVecs(H, NULL, &xs); CHKERRQ(ierr);
  ierr = MatCreateVecs(H, NULL, &ys); CHKERRQ(ierr);

  // Eigensolver
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, H, S); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  // Solve 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  ierr = EPSGetIterationNumber(eps, &its); CHKERRQ(ierr);
  ierr = EPSGetType(eps, &type); CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps, &nev, NULL, NULL);
  ierr = EPSGetTolerances(eps, &tol, &maxit); CHKERRQ(ierr);
  
  ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
  for(i = 0; i < nconv; i++) {
    ierr = EPSGetEigenpair(eps, i, &eig, &im_eig, xs, ys); CHKERRQ(ierr);
    ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%f (%g)", eig, error); CHKERRQ(ierr);
  }

  // Finalization
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&H); CHKERRQ(ierr);
  ierr = MatDestroy(&S); CHKERRQ(ierr);
  ierr = VecDestroy(&xs); CHKERRQ(ierr);
  ierr = VecDestroy(&ys); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}
