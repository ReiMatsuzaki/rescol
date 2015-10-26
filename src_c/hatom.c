#include "mat.h"
#include "fem_inf.h"

static char help[] = "solve hydrogen atom eigenvalue problem";

/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_mpd : maximum dim of projected problem
  
  -fem_type
  -bss_order or -dvr_nq
  -bps_type
  -bps_num_zs
  -bps_zmax
*/

int main(int argc, char **args) {

  PetscErrorCode ierr;
  FEMInf fem;
  MPI_Comm comm = PETSC_COMM_SELF;
  char target_dir[100]  = ".";
  char guess_type[10] = "none";
  PetscInt L = 0;
  time_t t0, t1, t2;
  
  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", &t0);
  PetscOptionsBegin(comm, "", "solve h atom eigen problem", "none");
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  PetscOptionsGetInt(NULL, "-L", &L, NULL);
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);
  PetscOptionsEnd();

  // Matrix
  Mat H, V, S;
  PrintTimeStamp(comm, "D2Mat", NULL);
  FEMInfSetD2R1Mat(fem, &H); MatScale(H, -0.5);
  PrintTimeStamp(comm, "VMat", NULL);
  FEMInfSetENR1Mat(fem, 0, 0.0, &V); MatAXPY(H, -1.0, V, SUBSET_NONZERO_PATTERN);

  PetscBool s_is_id; FEMInfGetOverlapIsId(fem, &s_is_id);
  if(s_is_id) 
    S = NULL;
  else {
    PrintTimeStamp(comm, "SMat", NULL);
    FEMInfSetSR1Mat(fem, &S);
  }

  // Solve
  EPS eps; 
  PrintTimeStamp(comm, "EPS", &t1);
  EPSCreate(comm, &eps);
  EPSSetOperators(eps, H, S);
  if(s_is_id)
    EPSSetProblemType(eps, EPS_HEP);
  else
    EPSSetProblemType(eps, EPS_GHEP);
  if(strcmp(guess_type, "heig") == 0) {
    PrintTimeStamp(comm, "guess", NULL);
    Vec guess[1]; 
    ierr = FEMInfGuessHEig(fem, 1, 0, &guess[0]); CHKERRQ(ierr);
    ierr = EPSSetInitialSpace(eps, 1, guess); CHKERRQ(ierr);
  }
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6);
  ierr = EPSSolve(eps);  CHKERRQ(ierr);

  // Output
  PrintTimeStamp(comm, "Output", &t2);
  FEMInfFPrintf(fem, stdout, 0);
  PetscPrintf(comm, "\n==== Time ====\n");
  PetscPrintf(comm, "t(mat): %f\n", t1-t0);
  PetscPrintf(comm, "t(diag): %f\n", t2-t1);

  PetscPrintf(comm, "\n==== EigenValues ====\n");
  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, NULL, NULL, NULL);
    PetscPrintf(comm, "eig%i: %f (%f)\n", i, kr, -1.0/((i+1)*(i+1)*2));
  }

  // Finalize
  MatDestroy(&H); MatDestroy(&S); MatDestroy(&V);
  EPSDestroy(&eps);
  FEMInfDestroy(&fem);
  return 0;
}
