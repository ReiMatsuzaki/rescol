#include <time.h>
#include <petscmat.h>
#include "mat.h"

static char help[] = "Read h_mat and s_mat written in COO format and solve eigen value problem";

/*
  important options in PETSc/SLEPc 
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_mpd : maximum dim of projected problem
  
  -eps_real : somputing nearest eigenvalues in real part
  -eps_target : target eigen value
*/

int AddWriteTime(const char* path_file, double t_read, double t_compute) {
  FILE* fp = NULL;
  fp = fopen(path_file, "a");
  PetscFPrintf(PETSC_COMM_WORLD, fp, "t_read: %f\n", t_read);
  PetscFPrintf(PETSC_COMM_WORLD, fp, "t_diag: %f\n", t_compute);
  fclose(fp);
  return 0;
}

int main (int argc, char **args)
{
  PetscErrorCode ierr;
  char h_path[256] = "hmat.dat";
  char s_path[256] = "smat.dat";
  Mat H, S;
  EPS eps;
  time_t t0, t1, t2, timer;
  struct tm *time0;

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
  t0 = time(&timer); 
  time0 = localtime(&timer);
  PetscPrintf(PETSC_COMM_WORLD, "t0 = %d h : %d m : %d s\n", time0->tm_hour, time0->tm_min, time0->tm_sec);

  MatCreateFromCOOFormatFile(h_path, &H);
  MatCreateFromCOOFormatFile(s_path, &S);
  t1 = time(&timer);
  time0 = localtime(&timer);
  PetscPrintf(PETSC_COMM_WORLD, "t1 = %d h : %d m : %d s\n", time0->tm_hour, time0->tm_min, time0->tm_sec);

  // Eigensolver
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, H, S); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_GHEP); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  // Solve 
  ierr = EPSSolve(eps); CHKERRQ(ierr);
  t2 = time(&timer);
  time0 = localtime(&timer);
  PetscPrintf(PETSC_COMM_WORLD, "t2 = %d h : %d m : %d s\n", time0->tm_hour, time0->tm_min, time0->tm_sec);

  // write to file
  EPSWriteToFile(eps, "diag_detail.dat", "diag_eigvals.dat", NULL);
  AddWriteTime("diag_detail.dat", 
	       (double)difftime(t1, t0),
	       (double)difftime(t2, t1));

  // Finalization
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&H); CHKERRQ(ierr);
  ierr = MatDestroy(&S); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}
