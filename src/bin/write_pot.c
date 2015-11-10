#include <rescol/pot.h>
#include <rescol/viewerfunc.h>

static char help[] = "write potential curve";

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;

  POT pot;
  ViewerFunc  viewer; ViewerFuncCreate(&viewer, comm);
  int J = 0;
  PetscReal mu = 1.0;

  ierr = PetscInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  ierr = PetscOptionsBegin(comm, "", "write_pot options", "none");

  ierr = POTCreateFromOptions(&pot, comm); CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(viewer); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-J", &J, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, "-mu", &mu, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if(!ViewerFuncIsActive(viewer))
    SETERRQ(comm, 1, "ViewerFunc is bad state");

  PetscPrintf(comm, "J=%d\n", J);
  PetscPrintf(comm, "mu=%f\n", mu);
  POTView(pot); 
  ViewerFuncView(viewer, PETSC_VIEWER_STDOUT_SELF);

  if(!ViewerFuncIsActive(viewer))
    POTViewFunc(pot, viewer);

  ierr = POTDestroy(&pot); CHKERRQ(ierr);
  ierr = ViewerFuncDestroy(&viewer); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

