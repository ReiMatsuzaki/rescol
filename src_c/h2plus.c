#include "mat.h"
#include "fem_inf.h"

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
  

  -target_dir : working directory. input and output files are located in this.
*/

PetscErrorCode MatCreateDirFile(const char* dn, const char* fn, Mat *M) {
  PetscErrorCode ierr;
  char path[100];
  sprintf(path, "%s/%s", dn, fn);
  ierr = MatCreateFromCOOFormatFile(path, M); CHKERRQ(ierr);  
  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  FEMInf fem;
  MPI_Comm comm = PETSC_COMM_SELF;
  char target_dir[100] = ".";
  char guess_type[10] = "none";
  double bond_length = 2.0;
  int qmax = 10;

  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &qmax, NULL);
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);
  PetscOptionsEnd();  

  // other conf
  PetscBool s_is_id; FEMInfGetOverlapIsId(fem, &s_is_id);

  // Matrix
  Mat H;

  // d2 term
  Mat d2_r1, s_y1;
  PrintTimeStamp(comm, "D2", NULL);
  ierr = FEMInfSetD2R1Mat(fem, &d2_r1); CHKERRQ(ierr);
  ierr = MatCreateDirFile(target_dir, "s_y1mat.dat", &s_y1); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(d2_r1, s_y1, comm, &H); CHKERRQ(ierr);
  ierr = MatScale(H, -0.5); CHKERRQ(ierr);
  ierr = MatDestroy(&d2_r1); CHKERRQ(ierr);
  
  // Overlap
  Mat S;
  if(s_is_id)
    S = NULL;
  else {
    Mat s_r1; 
    PrintTimeStamp(comm, "S", NULL);    
    ierr = FEMInfSetSR1Mat(fem, &s_r1); CHKERRQ(ierr);
    ierr = MatSetSynthesizeFast(s_r1, s_y1, comm, &S);  CHKERRQ(ierr);
    ierr = MatDestroy(&s_r1); CHKERRQ(ierr);
  }
  ierr = MatDestroy(&s_y1);  CHKERRQ(ierr);

  // angular term
  Mat L, l_r1, l_y1;
  PrintTimeStamp(comm, "L", NULL); CHKERRQ(ierr);
  ierr = FEMInfSetR2invR1Mat(fem, &l_r1); CHKERRQ(ierr);
  ierr = MatCreateDirFile(target_dir, "l_y1mat.dat", &l_y1); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(l_r1, l_y1, comm, &L); CHKERRQ(ierr);
  ierr = MatDestroy(&l_r1);  CHKERRQ(ierr);
  ierr = MatDestroy(&l_y1); CHKERRQ(ierr);
  ierr = MatAXPY(H, 0.5, L, SUBSET_NONZERO_PATTERN);

  // Nuclear-Electron term
  PrintTimeStamp(comm, "N-E", NULL);
  PetscReal a = bond_length / 2.0;
  for(int q = 0; q < qmax; q++) {
    char pq_path[100]; sprintf(pq_path, "%s/p%d_y1mat.dat", target_dir, q);
    FILE* f = fopen(pq_path, "r");
    if(f != NULL) {
      Mat q_y1, q_r1, V;
      ierr = FEMInfSetENR1Mat(fem, q, a, &q_r1); CHKERRQ(ierr);
      ierr = MatCreateFromCOOFormatFileHandler(f, &q_y1); CHKERRQ(ierr);
      ierr = MatSetSynthesizeFast(q_r1, q_y1, comm,  &V); CHKERRQ(ierr);
      ierr = MatAXPY(H, -2.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = MatDestroy(&q_y1); CHKERRQ(ierr);
      ierr = MatDestroy(&q_r1); CHKERRQ(ierr);
      ierr = MatDestroy(&V); CHKERRQ(ierr);
      fclose(f);
    }
  }

  // Solve
  PrintTimeStamp(comm, "EPS", NULL);
  EPS eps; 
  EPSCreate(comm, &eps);
  EPSSetOperators(eps, H, S);
  if(s_is_id)
    EPSSetProblemType(eps, EPS_HEP);
  else
    EPSSetProblemType(eps, EPS_GHEP);
  EPSSetTarget(eps, -4.0);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSolve(eps);  
  
  // Output
  PrintTimeStamp(comm, "Output", NULL);
  PetscPrintf(PETSC_COMM_SELF, "START H atom\n");
  FEMInfFPrintf(fem, stdout, 0);

  int nconv;
  PetscScalar kr, ki;
  Vec xr, xi;
  MatCreateVecs(H, NULL, &xr); MatCreateVecs(H, NULL, &xi);
  EPSGetConverged(eps, &nconv);
  PetscPrintf(PETSC_COMM_SELF, "\n==== Eigenvalues ====\n");
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
  }  
  return 0;
}

