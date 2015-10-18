#include "mat.h"
#include "bspline.h"

static char help[] = "create initial guess for He atom from hydrogen eigen function";
/*
  ======= EPS =========
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  ======= BSS ==========
  -bss_order :
  -bss_rmax :
  -bss_knots_num :

  ======= other ==========
  -l1 :
  -l2 :
  -in_dir
  -out_dir
*/

PetscErrorCode CalcHVec(BSS bss, PetscReal z, PetscInt L, MPI_Comm comm, Vec *x, PetscScalar *e) {
  PetscErrorCode ierr;
  
  Mat H;
  BSSSetD2R1Mat(bss, comm, &H); MatScale(H, -0.5);
  if(L!=0) {
    Mat R2inv;
    ierr = BSSSetR2invR1Mat(bss, comm, &R2inv); CHKERRQ(ierr);
    MatAXPY(H, L*(L+1)*0.5, R2inv, DIFFERENT_NONZERO_PATTERN); 
    MatDestroy(&R2inv);
  }
  Mat V;
  ierr = BSSSetENR1Mat(bss, 0, 0.0, comm, &V); CHKERRQ(ierr);
  MatAXPY(H, -z, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V);

  Mat S;
  ierr = BSSSetSR1Mat(bss, comm, &S); CHKERRQ(ierr);
  
  EPS eps;
  EPSCreate(comm, &eps); EPSSetOperators(eps, H, S);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6); EPSSetProblemType(eps, EPS_GHEP);
  EPSSetFromOptions(eps); 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  Vec xi;
  PetscScalar ei;
  PetscInt n;
  ierr = MatCreateVecs(S, NULL, x); CHKERRQ(ierr);
  ierr = MatCreateVecs(S, NULL, &xi); CHKERRQ(ierr);
  EPSGetConverged(eps, &n);
  if (n == 0) {
    SETERRQ(PETSC_COMM_WORLD, 1, "Calculation Failed");
  }
  EPSGetEigenpair(eps, 0, e, &ei, *x, xi);
  
  VecDestroy(&xi);
  EPSDestroy(&eps);
  MatDestroy(&S); MatDestroy(&H); 

  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  BSS bss;
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscReal z = 2.0;
  PetscInt L1 = 0;
  PetscInt L2 = 0;
  char in_dir[100] = ".";
  char out_dir[100] = ".";
  
  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PetscPrintf(comm, "\nhe_guess start\n");
  time_t t0; PrintTimeStamp(comm, "Init", &t0 );
  PetscOptionsBegin(comm, "", "he_guess.c options", "none");
  PetscOptionsGetString(NULL, "-in_dir", in_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-out_dir", out_dir, 100, NULL);
  PetscOptionsGetInt(NULL, "-L1", &L1, NULL);
  PetscOptionsGetInt(NULL, "-L2", &L2, NULL);
  PetscOptionsGetInt(NULL, "-L2", &L2, NULL);
  PetscOptionsGetReal(NULL, "-z", &z, NULL);
  ierr = BSSCreateFromOptions(&bss, comm); CHKERRQ(ierr);
  PetscOptionsEnd();  

  // Calculate Hydrogen
  PrintTimeStamp(comm, "H atom", NULL);
  Vec x1, x2;
  PetscScalar e1, e2;
  ierr = CalcHVec(bss, z, L1, comm, &x1, &e1); CHKERRQ(ierr);
  ierr = CalcHVec(bss, z, L2, comm, &x2, &e2); CHKERRQ(ierr);

  // y2
  Vec y2;
  char path[100]; sprintf(path, "%s/guess_y2vec.dat", in_dir);
  ierr = VecCreateFromFile(path, comm, &y2);

  // r2y2
  Vec r2, r2y2;
  ierr = VecSetSynthesize(x1, x2, 1.0, comm, &r2); CHKERRQ(ierr);
  ierr = VecSetSynthesize(r2, y2, 1.0, comm, &r2y2); CHKERRQ(ierr);
  
  // output
  PrintTimeStamp(comm, "Output", NULL);
  PetscViewer viewer;
  char out_path[100]; sprintf(out_path, "%s/guess.vec.dat", out_dir);
  PetscViewerBinaryOpen(comm, out_path, FILE_MODE_WRITE, &viewer);
  VecView(r2y2, viewer);

  PetscPrintf(comm, "z: %lf\n", z);
  PetscPrintf(comm, "L1: %d\n", L1);
  PetscPrintf(comm, "L2: %d\n", L2);
  PetscPrintf(comm, "in_dir: %s\n", in_dir);
  PetscPrintf(comm, "out_dir: %s\n", out_dir);
  PetscPrintf(comm, "eig1: %f\n", e1);
  PetscPrintf(comm, "eig2: %f\n", e2);
  BSSFPrintf(bss, comm, stdout, 0);
  
  PrintTimeStamp(comm, "Finalize", NULL);
  VecDestroy(&x1); VecDestroy(&x2); VecDestroy(&y2);
  VecDestroy(&r2); VecDestroy(&r2y2);
  PetscPrintf(comm, "he_guess end\n");

  SlepcFinalize();

  return 0;
}
