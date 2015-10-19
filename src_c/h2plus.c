#include "mat.h"
#include "bspline.h"

static char help[] = "solve H2^+ problem";
/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  -lmax :
  -gerade or -ungerade
  -bss_order : 
  -bss_rmax  : 
  -bss_knots_num :
  -bond_length : 

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
  BSS bss;
  MPI_Comm comm = PETSC_COMM_SELF;
  char target_dir[100] = ".";
  double bond_length = 2.0;
  int qmax = 10;

  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &qmax, NULL);
  ierr = BSSCreateFromOptions(&bss, comm);  CHKERRQ(ierr);
  PetscOptionsEnd();  

  // Matrix
  Mat H, S;

  // d2 term
  Mat d2_r1, s_y1;
  PrintTimeStamp(comm, "D2", NULL);
  ierr = BSSSetD2R1Mat(bss, comm, &d2_r1); CHKERRQ(ierr);
  ierr = MatCreateDirFile(target_dir, "s_y1mat.dat", &s_y1); CHKERRQ(ierr);
  ierr = MatInitSynthesize(d2_r1, s_y1, comm, &H); CHKERRQ(ierr);
  ierr = MatSynthesize(d2_r1, s_y1, -0.5, &H, ADD_VALUES); CHKERRQ(ierr);
  ierr = MatDestroy(&d2_r1); CHKERRQ(ierr);
  
  // Overlap
  Mat s_r1;
  PrintTimeStamp(comm, "S", NULL);
  ierr = BSSSetSR1Mat(bss, comm, &s_r1); CHKERRQ(ierr);
  ierr = MatInitSynthesize(s_r1, s_y1, comm, &S); CHKERRQ(ierr);
  ierr = MatSynthesize(s_r1, s_y1, 1.0, &S, INSERT_VALUES); CHKERRQ(ierr);
  ierr = MatDestroy(&s_r1); CHKERRQ(ierr);
  ierr = MatDestroy(&s_y1); CHKERRQ(ierr);

  // angular term
  Mat l_r1, l_y1;
  PrintTimeStamp(comm, "L", NULL); CHKERRQ(ierr);
  ierr = BSSSetR2invR1Mat(bss, comm, &l_r1); CHKERRQ(ierr);
  ierr = MatCreateDirFile(target_dir, "l_y1mat.dat", &l_y1); CHKERRQ(ierr);
  ierr = MatSynthesize(l_r1, l_y1, 0.5, &H, ADD_VALUES); CHKERRQ(ierr);
  ierr = MatDestroy(&l_r1);  CHKERRQ(ierr);
  ierr = MatDestroy(&l_y1); CHKERRQ(ierr);


  // Nuclear-Electron term
  PrintTimeStamp(comm, "N-E", NULL);
  PetscReal a = bond_length / 2.0;
  for(int q = 0; q < qmax; q++) {
    char pq_path[100]; sprintf(pq_path, "%s/p%d_y1mat.dat", target_dir, q);
    FILE* f = fopen(pq_path, "r");
    if(f != NULL) {
      Mat q_y1, q_r1;
      ierr = MatCreateFromCOOFormatFileHandler(f, &q_y1); CHKERRQ(ierr);
      BSSSetENR1Mat(bss, q, a, comm, &q_r1);
      MatSynthesize(q_r1, q_y1, -2.0, &H, ADD_VALUES);
      ierr = MatDestroy(&q_y1); CHKERRQ(ierr);
      ierr = MatDestroy(&q_r1); CHKERRQ(ierr);
      fclose(f);
    }
  }

  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);

  // Solve
  PrintTimeStamp(comm, "EPS", NULL);
  EPS eps; 
  EPSCreate(comm, &eps);
  EPSSetOperators(eps, H, S);
  EPSSetProblemType(eps, EPS_GHEP);
  EPSSetTarget(eps, -2.0);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSolve(eps);  
  
  // Output
    PrintTimeStamp(comm, "Output", NULL);
  PetscPrintf(PETSC_COMM_SELF, "START H atom\n");
  BSSFPrintf(bss, PETSC_COMM_SELF, stdout, 0);

  int nconv;
  PetscScalar kr, ki;
  Vec xr, xi;
  MatCreateVecs(H, NULL, &xr); MatCreateVecs(H, NULL, &xi);
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
  }  
  return 0;
}

