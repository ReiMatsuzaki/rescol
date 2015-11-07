#include <rescol/fem_inf.h>
#include <rescol/angmoment.h>

static char help[] = "solve H atom photoionization problem";

struct _p_HAtom {
  MPI_Comm comm;
  PetscInt L;
  FEMInf fem;
};
typedef struct _p_HAtom* HAtom;

PetscErrorCode CalcMat(FEMInf fem, int L, Mat *H, Mat *S) {

  PetscErrorCode ierr;
  char label[10]; sprintf(label, "L+%d", L);
  PrintTimeStamp(fem->comm, label, NULL);

  PetscBool s_is_id; FEMInfGetOverlapIsId(fem, &s_is_id);
  if(s_is_id)
    S = NULL;
  else {
    ierr = FEMInfSetSR1Mat(fem, S); CHKERRQ(ierr); CHKERRQ(ierr);
  }
  
  ierr = FEMInfSetD2R1Mat(fem, H); CHKERRQ(ierr);
  MatScale(*H, -0.5);

  if(L != 0) {
    Mat A;
    ierr = FEMInfSetR2invR1Mat(fem, &A); CHKERRQ(ierr);
    MatAXPY(*H, 0.5*L*(L+1), A, DIFFERENT_NONZERO_PATTERN);
  }

  Mat V;
  FEMInfSetENR1Mat(fem, 0, 0.0, &V);
  MatAXPY(*H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  return 0;
}
PetscErrorCode SolveInit(FEMInf fem, int L, PetscScalar *e0, Vec *x) {

  PetscErrorCode ierr;
  Mat H, S;
  ierr = CalcMat(fem, L, &H, &S); CHKERRQ(ierr);

  EPS eps;
  ierr = PrintTimeStamp(fem->comm, "EPS", NULL); CHKERRQ(ierr);
  ierr = EPSCreate(fem->comm, &eps); CHKERRQ(ierr);
  ierr = EPSSetTarget(eps, -0.6); CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);  CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, H, S); CHKERRQ(ierr);
  if(S == NULL) {
    ierr = EPSSetProblemType(eps, EPS_NHEP);  CHKERRQ(ierr);
  } else {
    ierr = EPSSetProblemType(eps, EPS_GNHEP);  CHKERRQ(ierr);
  }
  Vec x0[1]; MatCreateVecs(H, &x0[0], NULL); 
  int num; FEMInfGetSize(fem, &num);
  for(int i = 0; i < num; i++) {
    VecSetValue(x0[0], i, 0.5, INSERT_VALUES);
  }
  VecAssemblyBegin(x0[0]); VecAssemblyEnd(x0[0]);
  EPSSetInitialSpace(eps, 1, x0);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  ierr = EPSSolve(eps); CHKERRQ(ierr);
  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  if(nconv == 0) 
    SETERRQ(fem->comm, 1, "Failed to digonalize in init state\n");

  Vec x_ans;
  MatCreateVecs(H, &x_ans, NULL);
  EPSGetEigenpair(eps, 0, e0, NULL, x_ans, NULL);

  EPSDestroy(&eps);

  PetscScalar v[1]; PetscInt idx[1] = {1};
  VecGetValues(x_ans, 1, idx, v);
  PetscScalar scale_factor = v[0] / cabs(v[0]);
  VecScale( x_ans, 1.0/scale_factor);

  PetscScalar norm0;
  Vec Sx;  MatCreateVecs(S, &Sx, NULL); 
  MatMult(S, x_ans, Sx); VecDot(x_ans, Sx, &norm0);

  VecScale(x_ans, 1.0/sqrt(norm0));

  *x = x_ans;
  return 0;
}
PetscErrorCode SolveFinal(FEMInf fem, int L1, PetscScalar energy, 
			  Vec x0, Vec *x1, PetscScalar *alpha) {

  PetscErrorCode ierr;

  Mat S, L, D; 
  CalcMat(fem, L1, &L, &S);
  MatAXPY(L, -energy, S, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&S);

  PetscScalar mat_ele_cos;
  // <Y_10 | P_q(cos theta) | Y_00>
  mat_ele_cos = Y1ElePq(1, 1, 0,
			0,    0);
  if(getenv("SHOW_DEBUG"))
    printf("mat_ele_cos = %f\n", PetscRealPart(mat_ele_cos));

  POT dp_length; 
  ierr = POTPowerCreate(&dp_length, mat_ele_cos, 1.0); CHKERRQ(ierr);
  ierr = FEMInfSetPOTR1Mat(fem, dp_length, &D); CHKERRQ(ierr);

  Vec driv;
  MatCreateVecs(L, &driv, NULL);
  ierr = MatMult(D, x0, driv); CHKERRQ(ierr);
  ierr = MatDestroy(&D); CHKERRQ(ierr);

  KSP ksp;
  ierr = KSPCreate(fem->comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, L, L); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  ierr = MatCreateVecs(L, x1, NULL);

  ierr = KSPSolve(ksp, driv, *x1); CHKERRQ(ierr);
  ierr = VecDot(*x1, driv, alpha); CHKERRQ(ierr);

  // KSPView(ksp, PETSC_VIEWER_STDOUT_SELF);
  KSPDestroy(&ksp);
  MatDestroy(&L);
  VecDestroy(&driv);
  return 0;

}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  FEMInf fem;
  PetscReal w = 1.0;
  int L0 = 0;
  int L1 = 1;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  PrintTimeStamp(comm, "Init", NULL);
  PetscOptionsBegin(comm, "", "h2plus.c options", "none");
  ierr = PetscOptionsGetInt(NULL, "-L0", &L0, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-L1", &L1, NULL); CHKERRQ(ierr);  
  ierr = PetscOptionsGetReal(NULL, "-w", &w, NULL); CHKERRQ(ierr);  
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);  
  PetscOptionsEnd();

  int numx = 800;

  PetscReal *xs; PetscMalloc1(numx, &xs); 
  for(int i = 0; i < numx; i++)
    xs[i] = i * 0.1;
  FILE *fp; 
  
  Vec x0, x1;
  PetscScalar e0, alpha;
  ierr = SolveInit(fem, L0, &e0, &x0); CHKERRQ(ierr);
  FILE *fp0;
  ierr = PetscFOpen(comm, "tmp/h_pi_psi0.dat", "w", &fp0); CHKERRQ(ierr);
  ierr = FEMInfWritePsi(fem, xs, numx , x0, fp0); CHKERRQ(ierr);
  ierr = PetscFClose(comm, fp0); CHKERRQ(ierr);

  if(getenv("SHOW_DEBUG")) {
    printf("E0=%f\n", PetscRealPart(e0));
  }

  ierr = SolveFinal(fem, L1, e0+w, x0, &x1, &alpha); CHKERRQ(ierr);

  FEMInfView(fem);
  PetscPrintf(comm, "alpha: %f, %f\n", 
	      PetscRealPart(alpha), 
	      PetscImaginaryPart(alpha));

  ierr = PetscFOpen(comm, "tmp/h_pi_psi.dat", "w", &fp); CHKERRQ(ierr);
  ierr = FEMInfWritePsi(fem, xs, numx , x1, fp); CHKERRQ(ierr);
  ierr = PetscFClose(comm, fp); CHKERRQ(ierr);

  FEMInfDestroy(&fem);
  
  return 0;  
}
