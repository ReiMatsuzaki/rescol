#include "mat.h"
#include <petscmat.h>

static char help[] = "Read h_mat and s_mat written in COO format and solve eigen value problem";

/*
  important options in PETSc/SLEPc 
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_mpd : maximum dim of projected problem
  
  -eps_real : somputing nearest eigenvalues in real part
  -eps_target : target eigen value
*/

int main (int argc, char **args)
{
  PetscErrorCode ierr;
  char h_path[256] = "hmat.dat";
  char s_path[256] = "smat.dat";
  Mat H, S;
  EPS eps;

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

  // Eigensolver
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, H, S); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  // Solve 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  // write to file
  //  EPSWriteToFile(eps, "diag_detail.dat", "diag_eigvals.dat", "diag_eigvecs.dat");
    EPSWriteToFile(eps, "diag_detail.dat", "diag_eigvals.dat", NULL);

  // Finalization
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&H); CHKERRQ(ierr);
  ierr = MatDestroy(&S); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}
