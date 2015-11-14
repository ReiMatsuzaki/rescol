#include <rescol/viewerfunc.h>
#include <rescol/fem_inf.h>
#include <rescol/wavefunc.h>

static char help[] = "fitting function with various Finite Element Method";

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  ierr = PetscInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  KSP    ksp;         KSPCreate(comm, &ksp);
  WaveFunc wave_func; WaveFuncCreate(comm, &wave_func);
  FEMInf fem;         FEMInfCreate(comm, &fem); 
  PetscViewer viewer= PETSC_VIEWER_STDOUT_SELF;
  ViewerFunc view_func; ViewerFuncCreate(comm, &view_func);

  PetscViewerFormat format;

  ierr = PetscOptionsBegin(comm, "", "eig_one.c options", "none");
  ierr = WaveFuncSetFromOptions(wave_func); CHKERRQ(ierr);
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);
  //  ierr = PotSetFromOptions(pot);  CHKERRQ(ierr);  
  ierr = ViewerFuncSetFromOptions(view_func); CHKERRQ(ierr);
  ierr = PetscOptionsGetViewer(comm, NULL, "-viewer", &viewer, &format, NULL);
  PetscOptionsEnd();

  Vec c; FEMInfCreateVec(fem, 1, &c);
  ierr = FEMInfFit(fem, wave_func, ksp, c);        CHKERRQ(ierr);
  ierr = FEMInfViewFunc(fem, c, view_func); CHKERRQ(ierr);

  ierr = PFView(wave_func, viewer);      CHKERRQ(ierr);
  ierr = FEMInfView(fem, viewer); CHKERRQ(ierr);
  ierr = ViewerFuncView(view_func, viewer); CHKERRQ(ierr);  

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
