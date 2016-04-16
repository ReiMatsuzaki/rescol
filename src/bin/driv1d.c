#include "../../include/rescol/fem_inf.h"
#include "../../include/rescol/pot.h"
#include "../../include/rescol/viewerfunc.h"

static char help[] = "solve 1D driven type Shrodinger equation.";

/*
  Solve 1D driven type Schrodinger equation.
  (T + Vs + Vl - E) y = s;
*/

// ---- Range ----
typedef struct {
  MPI_Comm comm;
  PetscReal x0;
  PetscReal x1;
  PetscInt num;
} p_Range;
typedef p_Range* Range;
PetscErrorCode RangeCreate(MPI_Comm comm, Range *p_self) {
  Range self;
  PetscNew(&self);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode RangeDestroy(Range *p_self) {
  PetscErrorCode ierr;
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode RangeSetFromOptions(Range self, const char prefix[]) {
  PetscErrorCode ierr;
  char range_string[100];
  PetscBool find;

  //  char option_name[100] = "-";
  //  strcat(option_name, prefix);
  //  strcat(option_name, "-range");
  char option_name[100] = "-energy_range";
  ierr = PetscOptionsGetString(NULL, option_name, range_string, 100, &find);
  CHKERRQ(ierr);
  if(!find) {
    /*
    char msg[1000];
    printf("a\n");
    sprintf(msg, "option %s is not found.", option_name);
    printf("aa\n");
    */
    SETERRQ(self->comm, 1, "range option");
    printf("b\n");
  }

  // -- Count up number of ":"
  int num_sep = 0;
  for(int i = 0; i < strlen(range_string); i++) {
    if(range_string[i] == ':')
      num_sep++;
  }
  
  if(num_sep == 0) {
    self->x0 = atof(range_string);
    self->x1 = self->x0;
    self->num = 1;
    return 0;
  }
  
  if(num_sep == 2) {
    char *str_x0 = strtok(range_string, ":");
    char *str_x1 = strtok(NULL, ":");
    char *str_num = strtok(NULL, ":");
    self->x0 = atof(str_x0);
    self->x1 = atof(str_x1);
    self->num= atol(str_num);
    return 0;
  }

  SETERRQ(self->comm, 1, "Format invalid.");

}
PetscErrorCode RangeView(Range self, PetscViewer v) {
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);  

  if(iascii) {
    PetscViewerASCIIPrintf(v, "Range object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "x0:  %f\n", self->x0);
    PetscViewerASCIIPrintf(v, "x1:  %f\n", self->x1);
    PetscViewerASCIIPrintf(v, "num: %d\n", self->num);  
    PetscViewerASCIIPopTab(v);
  } else {
    SETERRQ(self->comm, 1, "only ascii is supported.");
  }
  return 0;
}
PetscReal RangeGetVal(Range self, int i ) {
  if(i < 0 || self->num <= i) {
    SETERRQ(self->comm, 1, "out of range");
  }
  if(i == 0)
    return self->x0;
  else 
    return self->x0 + i * (self->x1-self->x0)/(self->num-1);
}

MPI_Comm comm = MPI_COMM_SELF;

// ---- Input ----
char pot_type[10] = "single";
FEMInf      fem;
ViewerFunc  viewer;
Pot         pot_v0;
Pot         pot_v1;
Pot         driv;
Range energy_range;
KSP ksp;
Vec c0, c1, m;
Mat VL, VS, S;  

PetscErrorCode Driv1dCreate() {
  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Init", NULL);

  ierr = FEMInfCreate(comm, &fem); CHKERRQ(ierr);  
  ierr = ViewerFuncCreate(comm, &viewer); CHKERRQ(ierr);  
  ierr = PotCreate(comm, &pot_v1); CHKERRQ(ierr);  
  ierr = PotCreate(comm, &pot_v0); CHKERRQ(ierr);  
  ierr = PotCreate(comm, &driv); CHKERRQ(ierr);  
  ierr = RangeCreate(comm, &energy_range); CHKERRQ(ierr);
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode Driv1dSetFromOptions() {
  PetscErrorCode ierr;
  
  PetscBool find;

  PrintTimeStamp(comm, "Set", NULL);
  PetscOptionsBegin(comm, "", "driv1d.c options", "none");

  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);  
  ierr = ViewerFuncSetFromOptions(viewer); CHKERRQ(ierr);
  ierr = RangeSetFromOptions(energy_range, "energy"); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-pot_type", pot_type, 10, &find);
  CHKERRQ(ierr);

  if(strcmp(pot_type, "single") == 0) {
    ierr = PotSetFromOptions2(pot_v0, "v0");  CHKERRQ(ierr);
  } else if(strcmp(pot_type, "double") == 0) {
    ierr = PotSetFromOptions2(pot_v0, "v0");  CHKERRQ(ierr);
    ierr = PotSetFromOptions2(pot_v1, "v1"); CHKERRQ(ierr);
  } else {
    SETERRQ(comm, 1, "option value of -pot_type must be single or double");
  }

  ierr = PotSetFromOptions2(driv,  "driv"); CHKERRQ(ierr);

  PetscOptionsEnd();

  ierr = FEMInfCreateMat(fem, 1, &VL); CHKERRQ(ierr);
  ierr = FEMInfCreateMat(fem, 1, &VS); CHKERRQ(ierr);
  ierr = FEMInfCreateMat(fem, 1, &S);  CHKERRQ(ierr);

  ierr = VecCreate(comm, &c0);  CHKERRQ(ierr);
  ierr = VecSetType(c0, "seq"); CHKERRQ(ierr);

  ierr = VecCreate(comm, &c1); CHKERRQ(ierr);
  ierr = VecSetType(c1, "seq"); CHKERRQ(ierr);  
  
  return 0;
}
PetscErrorCode PrintIn() {

  PetscErrorCode ierr;

  // -- print input information --
  RangeView(energy_range, PETSC_VIEWER_STDOUT_SELF);
  if(strcmp(pot_type, "single") == 0) {
    ierr = PFView(pot_v0, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  } else {
    printf("V0 range potential:\n");
    ierr = PFView(pot_v0,  PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    printf("V1 range potential:\n");
    ierr = PFView(pot_v1, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  }
  return 0;
}
PetscErrorCode Driv1dCalc1(PetscReal energy) {

  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Calc1", NULL);

  // -- Discretise operator --
  // L1 = E - T - VL - VS
  // L0 = E - T - VL

  Mat L0;
  ierr = FEMInfCreateMat(fem, 1, &L0); CHKERRQ(ierr);  
  ierr = FEMInfD2R1Mat(fem, L0);       CHKERRQ(ierr);
  ierr = MatScale(L0, 0.5);            CHKERRQ(ierr);
  
  ierr = FEMInfPotR1Mat(fem, pot_v0, VL); CHKERRQ(ierr);
  ierr = MatAXPY(L0, -1.0, VL, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

  ierr = FEMInfSR1Mat(fem, S);        CHKERRQ(ierr);
  ierr = MatAXPY(L0, +energy,
		 S, DIFFERENT_NONZERO_PATTERN);  CHKERRQ(ierr);

  // -- Fitting --
  Vec c_fit;
  ierr = VecCreate(comm, &c_fit);  CHKERRQ(ierr); 
  ierr = VecSetType(c_fit, "seq"); CHKERRQ(ierr);
  ierr = FEMInfFit(fem, driv, ksp, c_fit); CHKERRQ(ierr);
  ierr = FEMInfCreateVec(fem, 1, &m); CHKERRQ(ierr);
  ierr = MatMult(S, c_fit, m);

  // -- Solve driven equation --
  ierr = KSPSetOperators(ksp, L0, L0); CHKERRQ(ierr);
  int n; FEMInfGetSize(fem, &n);
  //ierr = VecSetSizes(c0, PETSC_DECIDE, n); CHKERRQ(ierr);
  ierr = VecSetSizes(c0, n, n); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c0);  CHKERRQ(ierr);  

  ierr = VecDestroy(&c_fit); CHKERRQ(ierr);

  if(strcmp(pot_type, "double") == 0) {
    Mat L1;
    ierr = FEMInfCreateMat(fem, 1, &L1); CHKERRQ(ierr);  
    ierr = MatCopy(L0, L1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = FEMInfPotR1Mat(fem, pot_v1, VS); CHKERRQ(ierr);
    ierr = MatAXPY(L1, -1.0, VS, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, L1, L1); CHKERRQ(ierr);    
    //ierr = VecSetSizes(c1, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetSizes(c1, n, n); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, m, c1);CHKERRQ(ierr);  
    ierr = MatDestroy(&L1); CHKERRQ(ierr);
  }

  ierr = MatDestroy(&L0);    CHKERRQ(ierr);
  return 0;
}
PetscErrorCode Driv1dCalc2(PetscReal energy) {

  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Calc2", NULL);

  // -- Technical note
  // VecDot(x, y, val)  gives
  //  val = y^H x
  //
  // VecTDot(x, y, val) gives
  //  val = y^T x
  //

  if(strcmp(pot_type, "single") == 0) {
    PetscScalar alpha;
    ierr = VecTDot(m, c0, &alpha); CHKERRQ(ierr);
    printf("alpha : %f %f\n", creal(alpha), cimag(alpha));
  } else if(strcmp(pot_type, "double") == 0) {
    PetscScalar alpha;
    ierr = VecTDot(m, c1, &alpha); CHKERRQ(ierr);
    printf("alpha : %f %f\n", creal(alpha), cimag(alpha));


    PetscScalar plmx_psi0p;
    ierr = VecDot(c0, m, &plmx_psi0p); CHKERRQ(ierr);

    PetscScalar psi0p_s = conj(plmx_psi0p);
    printf("alpha0 : %f %f\n", creal(plmx_psi0p), cimag(plmx_psi0p));

    Vec v_c1; VecDuplicate(c1, &v_c1);
    ierr = MatMult(VS, c1, v_c1); CHKERRQ(ierr);
    PetscScalar psi0p_v_psip, psi0m_v_psip;
    ierr = VecDot(v_c1, c0, &psi0p_v_psip); CHKERRQ(ierr);
    ierr = VecTDot(v_c1,c0, &psi0m_v_psip); CHKERRQ(ierr);

    PetscReal k = sqrt(2.0 * energy);

    PetscScalar j_plmx = sqrt(-k * cimag(plmx_psi0p));

    PetscScalar impsi0p_v_psip = (psi0p_v_psip - psi0m_v_psip)/(2.0*I);
    PetscScalar amp = (cimag(psi0p_s) + impsi0p_v_psip) / j_plmx;
    printf("amplitude : %f, %f\n", creal(amp), cimag(amp));
    printf("phase : %f\n", carg(amp));
    printf("asb_amp : %f\n", cabs(amp));
    printf("abs_amp2_direct : %f\n", cabs(amp)*cabs(amp));
    printf("asb_amp_2_alpha : %f\n", -1.0/k * cimag(alpha));
  }

  return 0;
}
PetscErrorCode Driv1dWriteWaveFunc() {
  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Writing", NULL);
  if(strcmp(pot_type, "single")==0) {
    ierr = FEMInfViewFunc(fem, c0, viewer); CHKERRQ(ierr);
  } else {
    ierr = FEMInfViewFunc(fem, c1, viewer); CHKERRQ(ierr);
  }
  return 0;

}
PetscErrorCode Driv1dDestroy() {
  PrintTimeStamp(comm, "Finalize", NULL);
  FEMInfDestroy(&fem);
  ViewerFuncDestroy(&viewer);
  PFDestroy(&pot_v0);
  PFDestroy(&pot_v1);
  PFDestroy(&driv);
  VecDestroy(&c0); VecDestroy(&c1);
  MatDestroy(&S);
  KSPDestroy(&ksp); VecDestroy(&c1);
   return 0; 
 
}
int main(int argc, char **args) {

  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  ierr = Driv1dCreate();         CHKERRQ(ierr);
  ierr = Driv1dSetFromOptions(); CHKERRQ(ierr);
  ierr = PrintIn(); CHKERRQ(ierr);
  for(int i = 0; i < energy_range->num; i++) {
    PetscReal energy = RangeGetVal(energy_range, i);
    printf("energy : %f\n", energy);
    ierr = Driv1dCalc1(energy); CHKERRQ(ierr);
    ierr = Driv1dCalc2(energy); CHKERRQ(ierr);
  }
  ierr = Driv1dWriteWaveFunc(); CHKERRQ(ierr);
  ierr = Driv1dDestroy();       CHKERRQ(ierr);

  return 0;
}



