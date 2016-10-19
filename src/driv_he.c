#include <math.h>

#include "../../include/rescol/fem_inf.h"
#include "../../include/rescol/oce2.h"
#include "../../include/rescol/eeps.h"

static char help[] =  "Solve driven equation for two electron one center problem using numerical basis";

// ---- Input ----
PetscReal w;    // -- Photon energy --
PetscScalar e0; // -- Ionization energy --
Vec c0;         // -- Initial state vector --
Vec c1;         // -- Initial state vector --
FEMInf fem;     // -- Finite element method for initial and final state. --
OCE2 oce_i;      // -- Angular part of ionitial state --
OCE2 oce_f;      // -- Angular part of final state --

// ---- Intermediate ----
EEPS eps;       // -- Used to solve initial state --
KSP ksp;        // -- Used to solve final state --

// ---- Common ----
MPI_Comm comm = MPI_COMM_SELF;
PetscErrorCode ierr;

int Create() {

  PrintTimeStamp(comm, "Init", NULL);  
  
  ierr = OCE2Create(comm, &oce_i); CHKERRQ(ierr);
  ierr = OCE2Create(comm, &oce_f); CHKERRQ(ierr);
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = EEPSCreate(comm, &eps); CHKERRQ(ierr);
  
  return 0;

}
int Set() {
  FEMInf fem; ierr = FEMInfCreate(comm, &fem); CHKERRQ(ierr);

  // -- Set objects --  
  PrintTimeStamp(comm, "Sets", NULL);
  ierr = PetscOptionsGetReal(NULL, "-w", &w, NULL);
  ierr = PetscOptionsGetReal(NULL, "-e0", &e0, NULL);
  ierr = EEPSSetTarget(eps, -3.0); CHKERRQ(ierr);
  ierr = EEPSSetFromOptions(eps); CHKERRQ(ierr);
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);

  // -- FEMInf --
  FEMInf fem_i; FEMInfCreate(comm, &fem_i); 
  FEMInfSetFromOptions(fem_i);
  FEMInf fem_f; FEMInfDuplicate(fem_i, &fem_f); FEMInfCopy(fem_i, fem_f);

  // -- initial state --
  Y2s y2s_i;   ierr = Y2sCreate(comm, &y2s_i);   CHKERRQ(ierr);
  OCE2 oce2_i; ierr = OCE2Create(comm, &oce2_i); CHKERRQ(ierr);
  int Li = 0;
  ierr = PetscOptionsGetInt(NULL, "-Li", &Li, NULL); CHKERRQ(ierr);
  int Mi = 0;
  ierr = PetscOptionsGetInt(NULL, "-Mi", &Mi, NULL); CHKERRQ(ierr);
  int lmax_i = Li;
  ierr = PetscOptionsGetInt(NULL, "-lmax_i", &lmax_i, NULL); CHKERRQ(ierr);
  ierr = Y2sSetLM(y2s_i, Li, Mi, lmax_i); CHKERRQ(ierr);
  ierr = OCE2Set(oce_i, fem_i, y2s_i); CHKERRQ(ierr);

  // -- final state --
  Y2s y2s_f;   ierr = Y2sCreate(comm, &y2s_f);   CHKERRQ(ierr);
  OCE2 oce2_f; ierr = OCE2Create(comm, &oce2_f); CHKERRQ(ierr);
  int Lf = 0;
  ierr = PetscOptionsGetInt(NULL, "-Lf", &Lf, NULL); CHKERRQ(ierr);
  int Mf = 0;
  ierr = PetscOptionsGetInt(NULL, "-Mf", &Mf, NULL); CHKERRQ(ierr);
  int lmax_f = Lf;
  ierr = PetscOptionsGetInt(NULL, "-lmax_f", &lmax_f, NULL); CHKERRQ(ierr);
  ierr = Y2sSetLM(y2s_f, Lf, Mf, lmax_f); CHKERRQ(ierr);
  ierr = OCE2Set(oce_f, fem_f, y2s_f); CHKERRQ(ierr);

  return 0;
}
int CalcInit() {

  PrintTimeStamp(comm, "InitState", NULL);

  // -- matrix for initial state --
  Mat H;
  ierr = OCE2TMat(oce_i, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);
  ierr = OCE2PlusVneMat(oce_i, 0.0, 1.0, H); CHKERRQ(ierr);
  ierr = OCE2PlusVeeMat(oce_i, H); CHKERRQ(ierr);

  Mat S;
  PetscBool s_is_id;
  ierr = OCE2SMat(oce_i, MAT_INITIAL_MATRIX, &S, &s_is_id); CHKERRQ(ierr);

  if(s_is_id) {
    ierr = EEPSSetOperators(eps, H, NULL); CHKERRQ(ierr);
  } else {
    ierr = EEPSSetOperators(eps, H, S); CHKERRQ(ierr);
  }
  
  // -- Solve initial state --
  PrintTimeStamp(comm, "Solve", NULL);
  ierr = EEPSSolve(eps); CHKERRQ(ierr);
  ierr = MatCreateVecs(H, &c0, NULL); CHKERRQ(ierr);
  //  ierr = EPSGetEigenpair(eps->eps, 0, &e0, NULL, NULL, NULL); CHKERRQ(ierr);
  ierr = EEPSGetEigenvector(eps, 0, c0); CHKERRQ(ierr);
  PetscPrintf(comm, "E0 = (%f , %f)\n", PetscRealPart(e0), PetscImaginaryPart(e0));
  //  ierr = VecView(c0, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  return 0;
}
int CalcFinal() {
  
  PrintTimeStamp(comm, "FinalState", NULL);

  // -- matrix for final state --
  Mat L;
  ierr = OCE2TMat(oce_f, MAT_INITIAL_MATRIX, &L); CHKERRQ(ierr);
  ierr = OCE2PlusVneMat(oce_f, 0.0, 1.0, L); CHKERRQ(ierr);
  ierr = OCE2PlusVeeMat(oce_f, L); CHKERRQ(ierr);
  ierr = MatScale(L, -1.0); CHKERRQ(ierr);
  
  Mat S;
  PetscBool s_is_id;
  ierr = OCE2SMat(oce_f, MAT_INITIAL_MATRIX, &S, &s_is_id); CHKERRQ(ierr);
  if(s_is_id) {
    int nr, ny;
    MatDuplicate(L, MAT_DO_NOT_COPY_VALUES, &S);
    OCE2GetSizes(oce_f, &nr, &ny);
    for(int idx = 0; idx < nr*ny; idx++) {
      MatSetValue(S, idx, idx, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY); 
  }
  ierr = MatAXPY(L, w+e0, S, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  
  Mat Z;
  ierr = OCE2ZMat(oce_f, oce_i->y2s, MAT_INITIAL_MATRIX, &Z); CHKERRQ(ierr);
  Vec m; MatCreateVecs(L, &m, NULL);
  ierr = MatMult(Z, c0, m); CHKERRQ(ierr);

  PrintTimeStamp(comm, "solve", NULL);
  ierr = MatCreateVecs(L, &c1, NULL); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, L, L); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c1); CHKERRQ(ierr);

  PetscScalar alpha;
  VecTDot(c1, m, &alpha);
  PetscPrintf(comm, "alpha = (%f, %f)", PetscRealPart(alpha), PetscImaginaryPart(alpha));

  return 0;

}
int main(int argc, char **args) {

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  Create();
  Set();
  OCE2View(oce_i, PETSC_VIEWER_STDOUT_SELF);
  CalcInit();
  CalcFinal();
  

  // -- matrix for final state --
  //  PrintTimeStamp(comm, "Final", NULL);  
  //  Mat H1;
  // -- Print --
  //  PrintTimeStamp(comm, "Print", NULL);
  //PetscInt n; EPSGetConverged(eps, &n);
  
  

}
