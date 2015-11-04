#include <rescol/mat.h>
#include <rescol/fem_inf.h>
#include <rescol/angmoment.h>

static char help[] = "create initial guess for He atom from hydrogen eigen function";
/*
  ======= EPS =========
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  ======= FEM ==========
  -fem_type:
  -bps_zmax :
  -bps_num_zs:
  -bps_type: 

  ======= other ==========
  -l1 :
  -l2 :
  -in_dir
  -out_dir
  -guess_type
*/

PetscErrorCode CalcHVec(FEMInf this, PetscReal z, PetscInt L, MPI_Comm comm, Vec *x, PetscScalar *e) {

  PetscErrorCode ierr;
  PetscBool s_is_id; FEMInfGetOverlapIsId(this, &s_is_id);
  
  Mat H;
  FEMInfSetD2R1Mat(this, &H); MatScale(H, -0.5);
  if(L!=0) {
    Mat R2inv;
    ierr = FEMInfSetR2invR1Mat(this, &R2inv); CHKERRQ(ierr);
    MatAXPY(H, L*(L+1)*0.5, R2inv, DIFFERENT_NONZERO_PATTERN); 
    MatDestroy(&R2inv);
  }
  Mat V;
  ierr = FEMInfSetENR1Mat(this, 0, 0.0, &V); CHKERRQ(ierr);
  MatAXPY(H, -z, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V);

  EPS eps;
  EPSCreate(comm, &eps); 
  if(s_is_id) {
    EPSSetOperators(eps, H, NULL);
    EPSSetProblemType(eps, EPS_HEP);
  }
  else {
    Mat S;
    ierr = FEMInfSetSR1Mat(this, &S); CHKERRQ(ierr);
    EPSSetOperators(eps, H, S);
    EPSSetProblemType(eps, EPS_GHEP);
  }
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6); 
  EPSSetFromOptions(eps); 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  PetscInt n;
  ierr = MatCreateVecs(H, NULL, x); CHKERRQ(ierr);
  EPSGetConverged(eps, &n);
  if (n == 0) {
    SETERRQ(PETSC_COMM_WORLD, 1, "Calculation Failed");
  }
  EPSGetEigenpair(eps, 0, e, NULL, *x, NULL);
  
  EPSDestroy(&eps);
  MatDestroy(&H); 

  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  FEMInf fem;
  Y2s y2s;
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscReal z = 2.0;
  PetscInt L1 = 0;
  PetscInt L2 = 0;
  char in_dir[100] = ".";
  char out_dir[100] = ".";
  char guess_type[10] = "calc";
  
  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PetscPrintf(comm, "\nhe_guess start\n");
  time_t t0; PrintTimeStamp(comm, "Init", &t0 );
  PetscOptionsBegin(comm, "", "he_guess.c options", "none");
  PetscOptionsGetReal(NULL, "-z", &z, NULL);
  PetscOptionsGetInt(NULL, "-L1", &L1, NULL);
  PetscOptionsGetInt(NULL, "-L2", &L2, NULL);
  PetscOptionsGetString(NULL, "-in_dir", in_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-out_dir", out_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);
  ierr = Y2sCreateFromOptions(&y2s, comm); CHKERRQ(ierr);
  PetscOptionsEnd();  

  // Calculate Hydrogen
  Vec x1, x2;
  if(strcmp(guess_type, "calc") == 0) {
    PrintTimeStamp(comm, "calc H atom", NULL);
    PetscScalar e1, e2;
    ierr = CalcHVec(fem, z, L1, comm, &x1, &e1); CHKERRQ(ierr);
    ierr = CalcHVec(fem, z, L2, comm, &x2, &e2); CHKERRQ(ierr);
    PetscPrintf(comm, "eig1: %f\n", e1);
    PetscPrintf(comm, "eig2: %f\n", e2);    
  } else if (strcmp(guess_type, "eig") == 0) {
    PrintTimeStamp(comm, "get eig", NULL);
    ierr = FEMInfGuessHEig(fem, 1, L1, z, &x1);
    ierr = FEMInfGuessHEig(fem, 1, L2, z, &x2);
  } else {
    SETERRQ(comm, 1, "-guess_type <- {calc, eig}");
  }

  // y2
  Vec y2;
  PrintTimeStamp(comm, "Y2", NULL);
  ierr = Y2sSetGuessY2Vec(y2s, L1, L2, &y2); CHKERRQ(ierr);
  /*
  char path[100]; sprintf(path, "%s/guess_y2vec.dat", in_dir);
  ierr = VecCreateFromFile(path, comm, &y2);
  */

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

  FEMInfFPrintf(fem, stdout, 0);
  
  PrintTimeStamp(comm, "Finalize", NULL);
  VecDestroy(&x1); VecDestroy(&x2); VecDestroy(&y2);
  VecDestroy(&r2); VecDestroy(&r2y2);
  PetscPrintf(comm, "he_guess end\n");

  SlepcFinalize();

  return 0;
}
