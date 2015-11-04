#include "fem_inf.h"

static char help[] = "solve H atom photoionization problem";

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  FEMInf fem;
  char in_dir[10] = ".";
  int L0 = 0;
  int L1 = 1;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  PetscOptionsGetString(NULL, "-in_dir", in_dir, 10, NULL);
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);  
  PetscOptionsEnd();

  PrintTimeStamp(this->comm, "Read", NULL);
  Vec x0;
  ierr = VecCreate(comm, &x0); CHKERRQ(ierr);
  PetscViewer viewer;
  char path[20]; sprintf(path, "%s/x0.vec.dat", this->in_dir);
  ierr = PetscViewerBinaryOpen(comm, path, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  ierr = VecLoad(x0, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  
  // initial state
  PrintTimeStamp(comm, "InitState", NULL);  
  Mat M, L, V, S;
  ierr = FEMInfSetD2R1Mat(fem, &M); MatScale(M, -0.5);
  ierr = FEMInfSetR2invR1Mat(fem, &L); 
  MatAXPY(M, +0.5*L1*(L1+1), L, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&L);
  ierr = FEMInfSetENR1Mat(fem, 0, 0.0, &V); 
  MatAXPY(M, -1.0, V, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&V);
  ierr = FEMInfSetSR1Mat(fem, &S); CHKERRQ(ierr);
  MatAXPY(M, -0.5, S, DIFFERENT_NONZERO_PATTERN);

  
  
  

  return 0;  
}
