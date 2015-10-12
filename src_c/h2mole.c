#include <time.h>
#include "mat.h"
#include "bspline.h"


static char help[] = "solve H2 mole eigen value problem using bspline basis";
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

PetscErrorCode MatSetDirFile(const char* dn, const char* fn, Mat *M) {
  PetscErrorCode ierr;
  char path[100];
  sprintf(path, "%s/%s", dn, fn);
  ierr = MatCreateFromCOOFormatFile(path, M); CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode PrintTimeStamp(MPI_Comm comm, const char* label) {
  time_t t; time(&t);
  PetscPrintf(comm, "[%10s] %s", label, ctime(&t));
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
  PrintTimeStamp(PETSC_COMM_SELF, "Init");
  PetscOptionsBegin(comm, "", "h2mole.c options", "none");
  PetscOptionsGetString(NULL, "-target_dir", target_dir, 100, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &qmax, NULL);
  ierr = BSSCreateFromOptions(&bss); CHKERRQ(ierr);
  PetscOptionsEnd();

  Mat s_r1, s_y2, S;
  PrintTimeStamp(PETSC_COMM_SELF, "SMat");
  ierr = BSSSetSR1Mat(bss, comm,                &s_r1); CHKERRQ(ierr);
  ierr = MatSetDirFile(target_dir, "s_y2mat.dat", &s_y2);   CHKERRQ(ierr);
  ierr = MatSetSynthesize3(s_r1, s_r1, s_y2, 1.0, comm, &S);

  Mat H;  
  ierr = MatInitSynthesize3(s_r1, s_r1, s_y2, comm, &H); CHKERRQ(ierr);
  
  // D2 
  Mat d2_r1;
  PrintTimeStamp(PETSC_COMM_SELF, "D2Mat");
  ierr = BSSSetD2R1Mat(bss, comm, &d2_r1); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = MatSynthesize3(s_r1, d2_r1, s_y2, -0.5, &H, ADD_VALUES); CHKERRQ(ierr);
  ierr = MatSynthesize3(d2_r1, s_r1, s_y2, -0.5, &H, ADD_VALUES); CHKERRQ(ierr);
  MatDestroy(&d2_r1);

  // L
  Mat l_r1, l_1_y2, l_2_y2;
  PrintTimeStamp(PETSC_COMM_SELF, "LMat");
  ierr = BSSSetR2invR1Mat(bss, comm, &l_r1); CHKERRQ(ierr);
  ierr = MatSetDirFile(target_dir, "l_1_y2mat.dat", &l_1_y2); CHKERRQ(ierr);
  ierr = MatSetDirFile(target_dir, "l_2_y2mat.dat", &l_2_y2); CHKERRQ(ierr);
  ierr = MatSynthesize3(s_r1, l_r1, l_1_y2, 0.5, &H, ADD_VALUES);
  ierr = MatSynthesize3(l_r1, s_r1, l_2_y2, 0.5, &H, ADD_VALUES);
  MatDestroy(&l_r1); MatDestroy(&l_1_y2); MatDestroy(&l_2_y2);

  // e-n
  PrintTimeStamp(PETSC_COMM_SELF, "ENMat");
  PetscScalar a = bond_length/2.0;
  for(int q = 0; q < qmax; q++) {
    char path1[100]; char path2[100]; 
    sprintf(path1, "%s/p%d_A1_y2mat.dat", target_dir, q);
    sprintf(path2, "%s/p%d_A2_y2mat.dat", target_dir, q);
    FILE* f1 = fopen(path1, "r"); FILE* f2 = fopen(path2, "r");
    if(f1 != NULL && f2 != NULL) {
      Mat pq_1_y2, pq_2_y2, q_r1;
      ierr = MatCreateFromCOOFormatFileHandler(f1, &pq_1_y2); CHKERRQ(ierr);
      ierr = MatCreateFromCOOFormatFileHandler(f2, &pq_2_y2); CHKERRQ(ierr);
      ierr = BSSSetENR1Mat(bss, q, a, comm, &q_r1);  CHKERRQ(ierr);
      ierr = MatSynthesize3(q_r1, s_r1, pq_1_y2,-2.0,&H,ADD_VALUES); CHKERRQ(ierr);
      ierr = MatSynthesize3(s_r1, q_r1, pq_2_y2,-2.0,&H,ADD_VALUES); CHKERRQ(ierr);
      MatDestroy(&pq_1_y2); MatDestroy(&pq_2_y2); MatDestroy(&q_r1);
    }
  }

  // e-e
  PrintTimeStamp(PETSC_COMM_SELF, "EEMat");
  for(int q = 0; q < qmax; q++) {
    char path[100]; sprintf(path, "%s/p%d_12_y2mat.dat", target_dir, q);
    FILE *f = fopen(path, "r"); 
    if(f != NULL) {
      Mat r2, y2; 
      ierr = MatCreateFromCOOFormatFileHandler(f, &y2); CHKERRQ(ierr);
      ierr = BSSSetEER2Mat(bss, q, comm, &r2); CHKERRQ(ierr);
      ierr = MatSynthesize(r2, y2, 1.0, &H, ADD_VALUES); CHKERRQ(ierr);
      MatDestroy(&r2); MatDestroy(&y2);
    }
  }

  // Finalize for matrix
  PrintTimeStamp(PETSC_COMM_SELF, "Assemble");
  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  MatDestroy(&s_r1); MatDestroy(&s_y2);

  // Solve
  EPS eps; 
  PrintTimeStamp(PETSC_COMM_SELF, "eps");
  EPSCreate(comm, &eps); EPSSetOperators(eps, H, S);
  EPSSetTarget(eps, -2.0); EPSSetProblemType(eps, EPS_GHEP); 
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  // Output
  PrintTimeStamp(PETSC_COMM_SELF, "Output");
  BSSFPrintf(bss, comm, stdout, 0);
  int nconv;
  PetscScalar kr, ki;
  Vec xr, xi;
  ierr = MatCreateVecs(H, NULL, &xr); CHKERRQ(ierr);
  ierr = MatCreateVecs(H, NULL, &xi); CHKERRQ(ierr);
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
  }  

  // Destroy
  VecDestroy(&xr); VecDestroy(&xi);
  EPSDestroy(&eps);
  MatDestroy(&H); MatDestroy(&S);

  return 0;
}
