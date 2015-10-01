#include "mat.hpp"

static char help[] = "Read h_mat and s_mat written in COO format and solve eigen value problem";

int main (int argc, char **args)
{

  PetscErrorCode ierr;
  
  // Initialization
  ierr = SlepcInitialize(&argc, &args, (char*)0, help);

  // Read from options
  char h_path[256];
  char s_path[256];
  PetscOptionsString("-h-mat", "path for H matrix data file", "", h_path, h_path, 256, NULL);
  PetscOptionsString("-s-mat", "path for S matrix data file", "", s_path, s_path, 256, NULL);

  Mat H, S;
  //  Vec xr, xi;
  //  PetscReal error, tol, re, im;
  PetscInt n = 10;
  PetscInt Istart, Iend, i;
  
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);

  for(i = Istart; i < Iend; i++) {
    if (i>0) { ierr = MatSetValue(A, i, i-1, -1.0, INSERT_VALUES); CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A, i, i+1, -1.0, INSERT_VALUES); CHKERRQ(ierr); }
    ierr = MatSetValue(A, i, i, 2.0, INSERT_VALUES); CHKERRQ(ierr); 
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  MatView(A, PETSC_VIEWER_STDOUT_SELF);

  ierr = MatDestroy(&A); CHKERRQ(ierr);

  //  ierr = MatCreateVecs(A, NULL, &xr); CHKERRQ(ierr);
  

  // Finalization
  //ierr = PetscFinalize();
  ierr = SlepcFinalize();
  return 0;

}
