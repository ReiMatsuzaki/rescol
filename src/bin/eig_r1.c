#include <rescol/fem_inf.h>

static char help[] = "solve one dimensional problem for resonance energy";

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  FEMInf fem;
  POT pot;
  
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(comm, "", "traj_one.c options", "none");
  FEMInfCreateFromOptions(&fem, comm);
  POTCreateFromOptions(&pot, comm);
  PetscOptionsEnd();

  Mat H;
  FEMInfSetD2R1Mat(fem, &H); MatScale(H, -0.5);

  Mat V;
  FEMInfSetPOTR1Mat(fem, pot, &V); MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);

  Mat S;
  FEMInfSEtSR1MatNullable(fem, &S); 

  EPS eps; EPSCreateForBoundState(&eps, comm, H, S, 0.1);

  EPSSolve(eps);

  int nconv;
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    PetscScalar k;
    EPSGetEigenpair(eps, i, &k, NULL, NULL, NULL);
#if defined(PETSC_USE_COMPLEX)
    PetscPrintf(comm, "%d %f %f\n", i, PetscRealPart(k), PetscImaginaryPart(k));
#else
    PetscPrintf(comm, "%d %f\n", i, PetscRealPart(k), PetscImaginaryPart(k));
#endif
  }
  
  EPSDestroy(&eps);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
