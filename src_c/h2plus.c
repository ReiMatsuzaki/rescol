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
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &qmax, NULL);
  ierr = BSSCreateFromOptions(&bss);  CHKERRQ(ierr);
  PetscOptionsEnd();  

  // Matrix
  Mat H, S, rmat, ymat;

  // d2 term
  BSSSetD2R1Mat(bss, comm, &rmat);
  ierr = MatCreateDirFile(target_dir, "s_y1mat.dat", &ymat); CHKERRQ(ierr);
  ierr = MatInitSynthesize(rmat, ymat, comm, &H); CHKERRQ(ierr);
  MatSynthesize(rmat, ymat, -0.5, &H, ADD_VALUES);
  MatDestroy(&rmat);
  
  // Overlap
  BSSSetSR1Mat(bss, comm, &rmat);
  MatInitSynthesize(rmat, ymat, comm, &S);
  MatSynthesize(rmat, ymat, 1.0, &S, INSERT_VALUES);
  MatDestroy(&rmat); MatDestroy(&ymat);

  // angular term
  BSSSetR2invR1Mat(bss, comm, &rmat);
  ierr = MatCreateDirFile(target_dir, "l_y1mat.dat", &ymat); CHKERRQ(ierr);
  ierr = MatSynthesize(rmat, ymat, 0.5, &H, ADD_VALUES); CHKERRQ(ierr);
  MatDestroy(&rmat); MatDestroy(&ymat);
 
  // Nuclear-Electron term
  PetscReal a = bond_length / 2.0;
  for(int q = 0; q < qmax; q++) {
    char pq_path[100]; sprintf(pq_path, "%s/p%d_y1mat.dat", target_dir, q);
    FILE* f = fopen(pq_path, "r");
    if(f != NULL) {
      PetscPrintf(PETSC_COMM_SELF, "(q=%d) included\n", q);
      ierr = MatCreateFromCOOFormatFileHandler(f, &ymat); CHKERRQ(ierr);
      BSSSetENR1Mat(bss, q, a, comm, &rmat);
      MatSynthesize(rmat, ymat, -2.0, &H, ADD_VALUES);
      MatDestroy(&rmat); MatDestroy(&ymat);
      fclose(f);
    }
  }
  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);

  // Solve
  EPS eps; 
  EPSCreate(comm, &eps);
  EPSSetOperators(eps, H, S);
  EPSSetProblemType(eps, EPS_GHEP);
  EPSSetTarget(eps, -2.0);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSolve(eps);  
  
  // Output
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

