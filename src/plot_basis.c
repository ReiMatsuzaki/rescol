/*
  Write basis data for plot. 2016/4/20
*/

#include "../../include/rescol/fem_inf.h"
#include "../../include/rescol/viewerfunc.h"


static char help[] = "plot basis data for plot.";

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  FEMInf fem;
  ViewerFunc viewer;
  int n;
  
  
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);  

  // ==== Initialize ====
  PrintTimeStamp(comm, "Init", NULL);
  ierr = FEMInfCreate(comm, &fem); CHKERRQ(ierr);
  ierr = ViewerFuncCreate(comm, &viewer); CHKERRQ(ierr);


  // ==== Set values ====
  PrintTimeStamp(comm, "Set", NULL);
  PetscOptionsBegin(comm, "", "plot_basis.c options", "none");
  
  ierr = FEMInfSetFromOptions(fem);        CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(viewer); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-n", &n, NULL); CHKERRQ(ierr);

  PetscOptionsEnd();


  // ==== Check input values ====
  int num; FEMInfGetSize(fem, &num);
  if(n < 0 || num <= n) {
    SETERRQ(comm, 1, "n msut be zero or positive integer and smaller than size of FEM");
  }
  
  
  // ==== Calc ====
  Vec c;
  ierr = VecCreateSeq(comm, num, &c); CHKERRQ(ierr);
  ierr = VecSet(c, 0.0); CHKERRQ(ierr);
  int indices[1] = {n};
  PetscScalar values[1] = {1.0};
  ierr = VecSetValues(c, 1, indices, values, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(c); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(c); CHKERRQ(ierr);

  FEMInfViewFunc(fem, c, viewer);

  // ==== Finalize ====

  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
