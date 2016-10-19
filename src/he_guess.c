#include <rescol/fem_inf.h>
#include <rescol/y2s.h>
#include <rescol/eeps.h>
#include <rescol/wavefunc.h>
#include <rescol/synthesize.h>

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

PetscErrorCode CalcHVec(FEMInf fem, PetscReal z, int L, Vec *x, PetscScalar *e) {

  PetscErrorCode ierr;
  PetscBool s_is_id; FEMInfGetOverlapIsId(fem, &s_is_id);
  
  Mat H; FEMInfCreateMat(fem, 1, &H); FEMInfD2R1Mat(fem, H); MatScale(H, -0.5);
  if(L!=0) {
    Mat R2inv; FEMInfCreateMat(fem, 1, &R2inv); FEMInfR2invR1Mat(fem, R2inv);
    MatAXPY(H, L*(L+1)*0.5, R2inv, DIFFERENT_NONZERO_PATTERN); 
    MatDestroy(&R2inv);
  }
  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V);
  MatAXPY(H, -z, V, DIFFERENT_NONZERO_PATTERN); 
  MatDestroy(&V);

  

  EEPS eps; EEPSCreate(fem->comm, &eps); 
  
  if(s_is_id) {
    EEPSSetOperators(eps, H, NULL);
  }
  else {
    Mat S; FEMInfCreateMat(fem, 1, &S); FEMInfSR1Mat(fem, S);
    EEPSSetOperators(eps, H, S);
  }
  EEPSSetTarget(eps, -0.6); 
  EEPSSetFromOptions(eps);
  
  ierr = EEPSSolve(eps); CHKERRQ(ierr);

  PetscInt n;
  ierr = MatCreateVecs(H, NULL, x); CHKERRQ(ierr);
  EPSGetConverged(eps->eps, &n);
  if (n == 0) {
    SETERRQ(PETSC_COMM_WORLD, 1, "Calculation Failed");
  }
  EPSGetEigenpair(eps->eps, 0, e, NULL, *x, NULL);
  
  EEPSDestroy(&eps);
  MatDestroy(&H); 

  return 0;
}
PetscErrorCode FitHVec(FEMInf fem, PetscReal z, int n, int L, Vec *x) {
  MPI_Comm comm =fem->comm;
  PetscErrorCode ierr;
  KSP ksp; KSPCreate(comm, &ksp);
  WaveFunc f; WaveFuncCreate(comm, &f); WaveFuncSetHEig(f, n, L, z);

  Vec xx;
  ierr = FEMInfCreateVec(fem, 1, &xx);
  ierr = FEMInfFit(fem, f, ksp, xx); CHKERRQ(ierr);
  *x = xx;

  KSPDestroy(&ksp);
  PFDestroy(&f);
  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  FEMInf fem;        
  Y2s y2s;           
  PetscReal z = 2.0;
  PetscInt L1 = 0;
  PetscInt L2 = 0;
  char guess_type[10] = "calc";
//  char out_path[256] = "";
  PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
  PetscViewer vec_viewer;
  PetscViewerFormat format;
  
  // Initialize
  PrintTimeStamp(comm, "Init", NULL );

  PetscOptionsGetReal(NULL, "-z", &z, NULL);
  PetscOptionsGetInt(NULL, "-L1", &L1, NULL);
  PetscOptionsGetInt(NULL, "-L2", &L2, NULL);
  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  PetscOptionsGetViewer(comm, NULL, "-vec_viewer", &vec_viewer, &format, NULL);
//  PetscOptionsGetString(NULL, "-out_path", out_path, 256, NULL);
  ierr = FEMInfCreate(comm, &fem);  CHKERRQ(ierr);
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);
  ierr = Y2sCreate(comm, &y2s);     CHKERRQ(ierr);
  ierr = Y2sSetFromOptions(y2s);    CHKERRQ(ierr);

  // Calculate Hydrogen
  Vec x1, x2;
  if(strcmp(guess_type, "calc") == 0) {
    PrintTimeStamp(comm, "calc H atom", NULL);
    PetscScalar e1, e2;
    ierr = CalcHVec(fem, z, L1, &x1, &e1); CHKERRQ(ierr);
    ierr = CalcHVec(fem, z, L2, &x2, &e2); CHKERRQ(ierr);
    PetscPrintf(comm, "eig1: %f\n", e1);
    PetscPrintf(comm, "eig2: %f\n", e2);    
  } else if (strcmp(guess_type, "fit") == 0) {
    PrintTimeStamp(comm, "fit H atom", NULL);
    ierr = FitHVec(fem, z, 1, L1, &x1); CHKERRQ(ierr);
    ierr = FitHVec(fem, z, 1, L2, &x2); CHKERRQ(ierr);
  } else {
    SETERRQ(comm, 1, "-guess_type <- {calc, fit}");
  }

  // y2
  Vec y2; Y2sCreateY2Vec(y2s, &y2);
  PrintTimeStamp(comm, "Y2", NULL);
  ierr = Y2sGuessY2Vec(y2s, L1, L2, y2); CHKERRQ(ierr);

  // r2y2
  Vec r2, r2y2;
  ierr = VecVecSynthesize(x1, x2, 1.0, MAT_INITIAL_MATRIX, &r2); CHKERRQ(ierr);
  ierr = VecVecSynthesize(r2, y2, 1.0, MAT_INITIAL_MATRIX, &r2y2); CHKERRQ(ierr);
  
  // output
  PrintTimeStamp(comm, "Output", NULL);
  VecView(r2y2, vec_viewer);
  PetscPrintf(comm, "z: %lf\n", z);
  PetscPrintf(comm, "L1: %d\n", L1);
  PetscPrintf(comm, "L2: %d\n", L2);
  FEMInfView(fem, viewer);

  PrintTimeStamp(comm, "Finalize", NULL);
  FEMInfDestroy(&fem);
  Y2sDestroy(&y2s);
  PetscViewerDestroy(&vec_viewer);
  VecDestroy(&x1); VecDestroy(&x2); VecDestroy(&y2);
  VecDestroy(&r2); VecDestroy(&r2y2);
  PetscPrintf(comm, "he_guess end\n");

  SlepcFinalize();

  return 0;
}
