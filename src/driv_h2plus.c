#include "../include/math.h"
#include "../include/range.h"
#include "../include/eeps.h"
#include "../include/oce1.h"

static char help[] = "solve driven equation for h2+ ion:\n      (E-T-V)psi=s\nby grid method.\n";

typedef struct {
  // -- Inputs --
  MPI_Comm comm;
  char path_in[100];       // binary Vec file for initial state
  char path_out[100];      // binary Vec file for final state
  PetscReal R;
  PetscReal  E0;           // energy of bound initial state
  int dip_q;               // dipole direction (-1,0,+1)
  ViewerFunc viewer_func;
  ViewerFunc viewer_func_v;
  Range      w_range;      // calculation w range  
  FEMInf fem;  // NOTE: fem is used in oce_0 and oce_1  
  Y1s y1s_0;   //  NOTE: y1s_0 is used in oce_0 and oce_1
  Y1s y1s_1;   //  NOTE: y1s_0 is used in oce_0 and oce_1
  OCE1 oce_0;  // initial state
  OCE1 oce_1;  // final state 1
  KSP ksp;     // solver
} p_DrivH2plus;
typedef p_DrivH2plus* DrivH2plus;

void PrintTitle() {

  char label[100];
  PetscBool set;
  
  // -- print compile time --
  printf("\n>>>> driv_h2plus program >>>>\n");
  printf("%s", help);
  PetscOptionsGetString(NULL, NULL, "-label",  label, 100, &set);
  if(set)
    printf("label: %s\n", label);
  else
    printf("label: nothing\n");
  printf("Compile_date: %s %s\n", __DATE__, __TIME__);
  
}
PetscErrorCode Create(MPI_Comm comm, DrivH2plus *p_self) {

  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Init", NULL);

  DrivH2plus self;
  ierr = PetscNew(&self); CHKERRQ(ierr);
  self->comm = comm;
  *p_self = self;

  self->R = 2.0;
  ierr = FEMInfCreate(comm, &self->fem); CHKERRQ(ierr);
  ierr = Y1sCreate(comm, &self->y1s_0); CHKERRQ(ierr);
  ierr = Y1sCreate(comm, &self->y1s_1); CHKERRQ(ierr);  
  ierr = OCE1Create(comm, &self->oce_0); CHKERRQ(ierr);
  ierr = OCE1Create(comm, &self->oce_1); CHKERRQ(ierr);
  strcpy(self->path_in, "");
  strcpy(self->path_out, "");
  ierr = ViewerFuncCreate(comm, &self->viewer_func); CHKERRQ(ierr);
  ierr = ViewerFuncCreate(comm, &self->viewer_func_v); CHKERRQ(ierr);  
  ierr = RangeCreate(comm, &self->w_range); CHKERRQ(ierr);
  ierr = KSPCreate(comm, &self->ksp); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode SetFromOptions(DrivH2plus self) {

  PetscErrorCode ierr;
  PetscBool set;
  PrintTimeStamp(self->comm, "Set", NULL);

  PetscOptionsBegin(self->comm, "", "driv1d.c options", "none");

  // -- system parameter --
  ierr = PetscOptionsGetReal(NULL, NULL, "-R", &self->R, &set); CHKERRQ(ierr);
  if(!set)
    SETERRQ(self->comm, 1, "-R is not set");
  ierr = PetscOptionsGetReal( NULL, NULL, "-E0", &self->E0, &set); CHKERRQ(ierr);
  if(!set)
    SETERRQ(self->comm, 1, "-E0 is not set");
  ierr = PetscOptionsGetInt(NULL,NULL, "-dip_q", &self->dip_q, &set); CHKERRQ(ierr);
  if(!set)
    SETERRQ(self->comm, 1, "-dip_q is not set");
  
  // -- FEMInf --
  ierr = FEMInfSetFromOptions(self->fem); CHKERRQ(ierr);

  // -- OCE1 --
  ierr = Y1sSetOptionsPrefix(self->y1s_0, "init_"); CHKERRQ(ierr);
  ierr = Y1sSetFromOptions(self->y1s_0); CHKERRQ(ierr);
  ierr = OCE1Set(self->oce_0, self->fem, self->y1s_0); CHKERRQ(ierr);
  
  ierr = Y1sSetOptionsPrefix(self->y1s_1, "final_"); CHKERRQ(ierr);
  ierr = Y1sSetFromOptions(self->y1s_1); CHKERRQ(ierr);
  ierr = OCE1Set(self->oce_1, self->fem, self->y1s_1); CHKERRQ(ierr);

  // -- path --  
  ierr = PetscOptionsGetString(NULL, NULL, "-in",  self->path_in, 100, &set);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-out", self->path_out, 100, &set);
  CHKERRQ(ierr);

  // -- set function viewer --
  ierr = ViewerFuncSetOptionsPrefix(self->viewer_func, "length_"); CHKERRQ(ierr);
  ierr = ViewerFuncSetOptionsPrefix(self->viewer_func_v, "velocity_"); CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(self->viewer_func); CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(self->viewer_func_v); CHKERRQ(ierr);
  
  // -- set w range --
  ierr = RangeSetFromOptions(self->w_range, "-w_range"); CHKERRQ(ierr);

  PetscOptionsEnd();
  
  // -- Inputs check --
  
  return 0;
  
}
PetscErrorCode PrintIn(DrivH2plus self) {

  PrintTimeStamp(self->comm, "PrintIn", NULL);
  PetscErrorCode ierr;
  PetscViewer v = PETSC_VIEWER_STDOUT_SELF;

  ierr = PetscPrintf(self->comm, "R: %f\n", self->R); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "E0: %f\n", self->E0); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "dip_q: %d\n", self->dip_q); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "in: %s\n", self->path_in); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "out: %s\n", self->path_out); CHKERRQ(ierr);
  if(ViewerFuncIsActive(self->viewer_func)) {
    ierr = PetscPrintf(self->comm, "viewer_func:", v); CHKERRQ(ierr);
    ierr = ViewerFuncView(self->viewer_func, v); CHKERRQ(ierr);
  }
  if(ViewerFuncIsActive(self->viewer_func_v)) {
    ierr = PetscPrintf(self->comm, "viewer_func_v:", v); CHKERRQ(ierr);
    ierr = ViewerFuncView(self->viewer_func_v, v); CHKERRQ(ierr);
  }  
  ierr = PetscPrintf(self->comm, "w_range:", v); CHKERRQ(ierr);
  ierr = RangeView(self->w_range, v); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "fem:"); CHKERRQ(ierr);
  ierr = FEMInfView(self->fem, v); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "y1s_0:", v); CHKERRQ(ierr);
  ierr = Y1sView(self->y1s_0, v); CHKERRQ(ierr);
  ierr = PetscPrintf(self->comm, "y1s_1:", v); CHKERRQ(ierr);
  ierr = Y1sView(self->y1s_1, v); CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode Calc(DrivH2plus self) {

  PetscErrorCode ierr;

  // -- Read initial state vector --
  PrintTimeStamp(self->comm, "Read c0", NULL);
  PetscViewer v_c0;
  ierr = PetscViewerBinaryOpen(self->comm, self->path_in, FILE_MODE_READ, &v_c0);
  CHKERRQ(ierr);
  Vec c0;
  ierr = OCE1CreateVec(self->oce_0, &c0); CHKERRQ(ierr);
  ierr = VecLoad(c0, v_c0); CHKERRQ(ierr);
    
  // -- Calculate final state Hamiltonian and Overlap--
  PrintTimeStamp(self->comm, "Mat", NULL);
  Mat H, S;
  ierr = OCE1TMat(self->oce_1, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);
  ierr = OCE1PlusVneMat(self->oce_1, self->R/2.0, 1.0, H); CHKERRQ(ierr);
  PetscBool s_is_id;
  ierr = OCE1SMat(self->oce_1, MAT_INITIAL_MATRIX, &S, &s_is_id); CHKERRQ(ierr);

  // -- Dipole matrix element --
  Mat Z;
  ierr = OCE1ZMat(self->oce_1, self->oce_0, self->dip_q,
		  MAT_INITIAL_MATRIX, &Z); CHKERRQ(ierr);
  Mat DZ;
  ierr = OCE1DZMat(self->oce_1, self->oce_0, self->dip_q,
		   MAT_INITIAL_MATRIX, &DZ); CHKERRQ(ierr);
  
  // -- driven term --
  Vec driv;     // length form
  ierr = OCE1CreateVec(self->oce_1, &driv); CHKERRQ(ierr);
  ierr = MatMult(Z, c0, driv); CHKERRQ(ierr);
  Vec driv_v;   // velocity form
  ierr = OCE1CreateVec(self->oce_1, &driv_v); CHKERRQ(ierr);
  ierr = MatMult(DZ, c0, driv_v); CHKERRQ(ierr);  

  // -- start loop over energies --
  PetscBool is_first = PETSC_TRUE;
  PetscReal w;
  while(RangeNext(self->w_range, &w)) {
    PetscReal ene = w + self->E0;
    PetscPrintf(self->comm, "(w, E) = (%f, %f)\n", w, ene);
    
    Mat L;
    ierr = MatConvert(H, MATSAME, MAT_INITIAL_MATRIX, &L); CHKERRQ(ierr);
    if(s_is_id) {
      ierr = MatShift(L, -ene); CHKERRQ(ierr);
    } else {
      ierr = MatAXPY(L, -ene, S, SUBSET_NONZERO_PATTERN); CHKERRQ(ierr);
    }
    ierr = MatScale(L, -1.0); CHKERRQ(ierr);
    ierr = KSPSetOperators(self->ksp, L, L); CHKERRQ(ierr);

    Vec c1;
    ierr = OCE1CreateVec(self->oce_1, &c1);  CHKERRQ(ierr);
    ierr = KSPSolve(self->ksp, driv, c1);    CHKERRQ(ierr);
    Vec c1_v;
    ierr = OCE1CreateVec(self->oce_1, &c1_v); CHKERRQ(ierr);
    ierr = KSPSolve(self->ksp, driv_v, c1_v); CHKERRQ(ierr);

    PetscScalar alpha;
    VecTDot(c1, driv, &alpha);
    alpha /= 3.0;
    printf("alpha() = %f, %f\n",
	   PetscRealPart(alpha),
	   PetscImaginaryPart(alpha));

    PetscScalar alpha_v;
    VecTDot(c1_v, driv_v, &alpha_v);
    alpha_v /= 3.0;
    printf("alpha(velocity) = %f, %f\n", PetscRealPart(alpha_v),
	   PetscImaginaryPart(alpha_v));

    PetscReal c_light = 137.035999139;
    PetscReal au2mb = 5.291772 * 5.291772;
      
    PetscReal cs;
    cs = -4.0 * M_PI * w / c_light * PetscImaginaryPart(alpha) * au2mb;    
    PetscReal cs_v;
    cs_v = -4.0 * M_PI / (w * c_light) * PetscImaginaryPart(alpha_v) * au2mb;
    printf("CrossSection(Mb) = %f, %f\n", cs, cs_v);

    if(is_first && ViewerFuncIsActive(self->viewer_func)) {
      ierr = OCE1ViewFunc(self->oce_1, c1, self->viewer_func); CHKERRQ(ierr);
      ierr = OCE1ViewFunc(self->oce_1, c1_v, self->viewer_func_v); CHKERRQ(ierr);
    }
    is_first = PETSC_FALSE;
    
    MatDestroy(&L); VecDestroy(&c1);
  }
  
  return 0;
}
PetscErrorCode Finalize(DrivH2plus self) {

  PrintTimeStamp(self->comm, "Finalize", NULL);
  ViewerFuncDestroy(&self->viewer_func);
  ViewerFuncDestroy(&self->viewer_func_v);
  RangeDestroy(&self->w_range);
  FEMInfDestroy(&self->fem);
  self->oce_0->fem = NULL;
  self->oce_1->fem = NULL;
  OCE1Destroy(&self->oce_0);
  OCE1Destroy(&self->oce_1);
  KSPDestroy(&self->ksp);
  return 0;
  
}
int main(int argc, char **args) {

  PetscErrorCode ierr;
  DrivH2plus driv;
  MPI_Comm comm = MPI_COMM_SELF;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  PrintTitle();
  ierr = Create(comm, &driv); CHKERRQ(ierr);
  ierr = SetFromOptions(driv); CHKERRQ(ierr);
  ierr = PrintIn(driv); CHKERRQ(ierr);
  ierr = Calc(driv); CHKERRQ(ierr);
  ierr = Finalize(driv); CHKERRQ(ierr);

  printf("<<<< driv_h2plus program <<<<\n\n");
  return 0;
  
}


