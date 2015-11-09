#include <rescol/oce1.h>

static char help[] = "solve H2^+ problem";
/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  -bond_length : 
  -lmax :
  -gerade or -ungerade

  -fem_type
  -fd_num  -fd_xmax
  -bss_order 
  -dvr_nq
  -bps_type -bps_num_zs -bps_zmax
*/

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  OCE1 oce1;
  //  char guess_type[10] = "none";
  double bond_length = 2.0;

  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  //  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  ierr = OCE1CreateFromOptions(&oce1, comm); CHKERRQ(ierr);
  PetscOptionsEnd();  

  // other conf
  PetscBool s_is_id; FEMInfGetOverlapIsId(oce1->fem, &s_is_id);

  // Matrix
  Mat H, S;
  PrintTimeStamp(comm, "Mat", NULL);
  ierr = OCE1SetTMat(oce1, &H); CHKERRQ(ierr);
  ierr = OCE1PlusVneMat(oce1, bond_length/2.0, 1.0, &H); CHKERRQ(ierr);
  ierr = OCE1SetSMatNullable(oce1, &S); CHKERRQ(ierr);

  // Solve
  EPS eps; 
  PrintTimeStamp(comm, "EPS", NULL);
  EPSCreateForBoundState(&eps, comm, H, S, -1.2);
  EPSSolve(eps);  
  
  // Output
  PrintTimeStamp(comm, "Output", NULL);
  OCE1View(oce1);

  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  PetscPrintf(PETSC_COMM_SELF, "\n==== Eigenvalues ====\n");
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, NULL, NULL, NULL);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
  }  
  return 0;
}

