#include <math.h>

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
// pot_type = single => single potential problem
//          = double => double potential problem and apply two potential formula.
// problem_type = driv    => solve driven equation.
//              = scatter => solve scattering problem.
char pot_type[10] = "single";
char problem_type[10] = "driv"; 
int         L;       // must
PetscReal   Z;       // charge
PetscBool   use_v0;
Pot         pot_v0;  // must
PetscBool   use_v1;
Pot         pot_v1;  // only when pot_type=double
Pot         driv;    // must

FEMInf      fem;     // must
Range energy_range;  // must
PetscBool   use_func_view; // if true, print out wave function
ViewerFunc  viewer;  // must

// ---- intermediate ----
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

  // -- set system type --
  ierr = PetscOptionsGetString(NULL, "-pot_type", pot_type, 10, &find);
  CHKERRQ(ierr);
  if(!find)
    SETERRQ(comm, 1, "-pot_type is not found"); 
  ierr = PetscOptionsGetString(NULL, "-problem_type", problem_type, 10, &find);
  CHKERRQ(ierr);
  if(!find)
    SETERRQ(comm, 1, "-problem_type is not found"); 

  // -- set L --
  ierr = PetscOptionsGetInt(NULL, "-L", &L, &find); CHKERRQ(ierr);
  if(!find)
    SETERRQ(comm, 1, "-L is not found"); 
  if(L < 0)
    SETERRQ(comm, 1, "L must be zero or positive integer"); 

  // -- set Z --
  ierr = PetscOptionsGetReal(NULL, "-Z", &Z, &find); CHKERRQ(ierr);
  if(!find)
    Z = 0.0;

  // -- set potential
  if(strcmp(pot_type, "single") == 0) {
    ierr = PotSetFromOptions2(pot_v0, "v0", &find);  CHKERRQ(ierr);
    use_v0 = find;
    use_v1 = PETSC_FALSE;
  } else if(strcmp(pot_type, "double") == 0) {
    ierr = PotSetFromOptions2(pot_v0, "v0", &find);  CHKERRQ(ierr);
    use_v0 = find;
    ierr = PotSetFromOptions2(pot_v1, "v1", &find); CHKERRQ(ierr);
    use_v1 = find;
  } else {
    SETERRQ(comm, 1, "option value of -pot_type must be single or double");
  }

  // -- set driven term --
  if(strcmp(problem_type, "driv") == 0) {
    ierr = PotSetFromOptions2(driv,  "driv", NULL); CHKERRQ(ierr);
  }
  else if(strcmp(problem_type, "scatter") == 0) {
    // driven term is determined when energy is determined because
    // in scattering calculation, driven term is energy dependent.
  }
  else {
    SETERRQ(comm, 1, "invalid option value for -problem_type");
  }
  // -- other --
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);  
  ierr = ViewerFuncSetFromOptions(viewer, &use_func_view); CHKERRQ(ierr);
  ierr = RangeSetFromOptions(energy_range, "energy"); CHKERRQ(ierr);

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

  // -- compile time --
  printf("Compile date: %s  %s\n", __DATE__, __TIME__);

  // -- print input information --
  printf("pot_type: %s\n", pot_type);
  printf("problem_type: %s\n", problem_type);
  printf("L: %d\n", L);
  printf("Z: %f\n", Z);
  printf("Energy range:\n");
  RangeView(energy_range, PETSC_VIEWER_STDOUT_SELF);
  FEMInfView(fem, PETSC_VIEWER_STDOUT_SELF);
  if(use_v0) {
    printf("V0 potential:\n");
    ierr = PFView(pot_v0,  PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  }
  if(use_v1) {
    printf("V1 potential:\n");
    ierr = PFView(pot_v1, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  }
  PFView(driv, PETSC_VIEWER_STDOUT_SELF);
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

  if(L != 0) {
    Pot LL;
    Mat LLMat;
    ierr = PotCreate(comm, &LL); CHKERRQ(ierr);
    ierr = PotSetPower(LL, L*(L+1)/(2.0), -2); CHKERRQ(ierr);
    ierr = FEMInfCreateMat(fem, 1, &LLMat); CHKERRQ(ierr);
    ierr = FEMInfPotR1Mat(fem, LL, LLMat); CHKERRQ(ierr);
    ierr = MatAXPY(L0, -1.0, LLMat, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatDestroy(&LLMat); CHKERRQ(ierr);
  }

  Pot PotVc; PotCreate(comm, &PotVc); PotSetPower(PotVc, -Z, -1.0);
  Mat Vc;    FEMInfCreateMat(fem, 1, &Vc); FEMInfPotR1Mat(fem, PotVc, Vc);
  ierr = MatAXPY(L0, -1.0, Vc, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatDestroy(&Vc); CHKERRQ(ierr);
  ierr = PFDestroy(&PotVc); CHKERRQ(ierr);

  if(use_v0) {
    ierr = FEMInfPotR1Mat(fem, pot_v0, VL); CHKERRQ(ierr);
    ierr = MatAXPY(L0, -1.0, VL, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  }
  ierr = FEMInfSR1Mat(fem, S);        CHKERRQ(ierr);
  ierr = MatAXPY(L0, +energy,
		 S, DIFFERENT_NONZERO_PATTERN);  CHKERRQ(ierr);

  // -- driven term --
  if(strcmp(problem_type, "scatter") == 0) {
    double k = sqrt(energy*2.0);
    Pot rbessel;
    ierr = PotCreate(comm, &rbessel); CHKERRQ(ierr);
    ierr = PotSetRBessel(rbessel, L, k); CHKERRQ(ierr);
    Pot pots[2];
    if(strcmp(pot_type, "single") == 0) {
      pots[0] = pot_v0; pots[1] = rbessel;
    } 
    else if(strcmp(pot_type, "double") == 0) {
      pots[0] = pot_v1; pots[1] = rbessel;
    }
    PotSetProduct(driv, 2, pots);
  }
  //printf("driven term\n");
  //  ierr = PFView(driv, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = FEMInfCreateVec(fem, 1, &m); CHKERRQ(ierr);
  ierr = FEMInfPotR1Vec(fem, driv, m);CHKERRQ(ierr);

  // -- Solve driven equation --
  ierr = KSPSetOperators(ksp, L0, L0); CHKERRQ(ierr);
  int n; FEMInfGetSize(fem, &n);
  ierr = VecSetSizes(c0, n, n); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, m, c0);  CHKERRQ(ierr);  

  if(strcmp(pot_type, "double") == 0) {
    Mat L1;
    ierr = FEMInfCreateMat(fem, 1, &L1); CHKERRQ(ierr);  
    ierr = MatCopy(L0, L1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = FEMInfPotR1Mat(fem, pot_v1, VS); CHKERRQ(ierr);
    ierr = MatAXPY(L1, -1.0, VS, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, L1, L1); CHKERRQ(ierr);    
    ierr = VecSetSizes(c1, n, n); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, m, c1);CHKERRQ(ierr);  
    ierr = MatDestroy(&L1); CHKERRQ(ierr);
  }

  ierr = MatDestroy(&L0);    CHKERRQ(ierr);
  return 0;
}
PetscErrorCode Driv1dCalc2_single_scatter(PetscReal energy) {

  /*
    A = 2mu/k * <jl(kr) | V(r) | Psi>
    . = 2mu/k * <jl(kr) | V(r) | Psi^sc>_R0
   */

  PetscErrorCode ierr;

  double k = sqrt(2.0*energy);
  Pot rbessel;
  ierr = PotCreate(comm, &rbessel);          CHKERRQ(ierr);
  ierr = PotSetRBessel(rbessel, L, k);       CHKERRQ(ierr);
  KSP ksp;
  ierr = KSPCreate(comm, &ksp);              CHKERRQ(ierr);
  int n; FEMInfGetSize(fem, &n);
  Vec jvec;
  ierr = VecCreate(comm, &jvec); CHKERRQ(ierr);
  ierr = VecSetType(jvec, VECSEQ); CHKERRQ(ierr);
  ierr = FEMInfFit(fem, rbessel, ksp, jvec); CHKERRQ(ierr);
  
  PetscScalar amp0;
    ierr = VecTDot(m, jvec, &amp0); CHKERRQ(ierr);
  PetscScalar amp1;
  ierr = VecTDot(m, c0,   &amp1); CHKERRQ(ierr);
  amp0 *= 2.0/k;
  amp1 *= 2.0/k;
  
  double phase = carg(amp0 + amp1);
  if(phase < 0.0)
    phase += M_PI;
  double cs = 4.0 * M_PI * sin(phase)*sin(phase)/(k*k);

  // -- surface form --
  PetscReal R0 = 40.3;
  PetscScalar y, dy;
  FEMInfPsiOne(fem, c0, R0, &y);
  FEMInfDerivPsiOne(fem, c0, R0, &dy);
  PetscScalar j, dj;
  j = sin(k*R0); dj = k*cos(k*R0);
  PetscScalar amp_surface = -1.0/k*(j*dy-y*dj);
  
  double phase_sur = carg( amp_surface);
  if(phase_sur < 0.0)
    phase_sur += M_PI;
  double cs_sur = 4.0*M_PI*sin(phase_sur)*sin(phase_sur)/(k*k);
  

  PFDestroy(&rbessel);
  VecDestroy(&jvec);

  printf("amp_int0  : %15.10f %15.10f\n", creal(amp0), cimag(amp0));
  printf("amp_int1  : %15.10f %15.10f\n", creal(amp1), cimag(amp1));
  printf("phase_int : %15.10f\n", phase);
  printf("amp_surfa : %15.10f %15.10f\n", creal(amp_surface), cimag(amp_surface));
  printf("phase_sur : %15.10f\n", phase_sur);

  printf("cross_sec_int : %15.10f\n", cs);
  printf("cross_sec_sur : %15.10f\n", cs_sur);

  
  return 0;
}
PetscErrorCode Driv1dCalc2_single_driv(PetscReal energy) {
  PetscErrorCode ierr;
  

  /*
  Vec s_m; VecDuplicate(m, &s_m); 
  ierr = MatMult(S, m, s_m); CHKERRQ(ierr);

  PetscScalar alpha;
  ierr = VecTDot(s_m, c0, &alpha); CHKERRQ(ierr);
  printf("c0Sm : %f %f\n", creal(alpha), cimag(alpha));
    */

  PetscScalar mc0;
  ierr = VecTDot(m, c0, &mc0); CHKERRQ(ierr);
  printf("c0.m  : %15.10f %15.10f\n", creal(mc0), cimag(mc0));
  return 0;
}
PetscErrorCode Driv1dCalc2_double_driv(PetscReal energy) {
  PetscErrorCode ierr;
  PetscScalar alpha;
  ierr = VecTDot(m, c1, &alpha); CHKERRQ(ierr);
  printf("alpha : %f %f\n", creal(alpha), cimag(alpha));


  PetscScalar plmx_psi0p;
  ierr = VecDot(c0, m, &plmx_psi0p); CHKERRQ(ierr);

  PetscScalar psi0p_s = conj(plmx_psi0p);
  printf("alpha0 : %15.10f %15.10f\n", creal(plmx_psi0p), cimag(plmx_psi0p));
  
  Vec v_c1; VecDuplicate(c1, &v_c1);
  ierr = MatMult(VS, c1, v_c1); CHKERRQ(ierr);
  PetscScalar psi0p_v_psip, psi0m_v_psip;
  ierr = VecDot(v_c1, c0, &psi0p_v_psip); CHKERRQ(ierr);
  ierr = VecTDot(v_c1,c0, &psi0m_v_psip); CHKERRQ(ierr);

  PetscReal k = sqrt(2.0 * energy);

  PetscScalar j_plmx = sqrt(-M_PI * cimag(plmx_psi0p));
  
  PetscScalar impsi0p_v_psip = (psi0p_v_psip - psi0m_v_psip)/(2.0*I);
  PetscScalar amp = (cimag(psi0p_s) + impsi0p_v_psip) / j_plmx;
  printf("amplitude : %15.10f, %15.10f\n", creal(amp), cimag(amp));
  printf("phase : %15.10f\n", carg(amp));
  printf("abs_amp : %15.10f\n", cabs(amp));
  printf("abs_amp2_direct : %15.10f\n", cabs(amp)*cabs(amp));
  printf("asb_amp_2_alpha : %15.10f\n", -1.0/k * cimag(alpha));
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

  if(strcmp(pot_type, "single") == 0 &&
     strcmp(problem_type, "scatter") == 0) {
    ierr = Driv1dCalc2_single_scatter(energy); CHKERRQ(ierr);
  }  

  if(strcmp(pot_type, "double") == 0 &&
     strcmp(problem_type, "scatter") == 0) {
    SETERRQ(comm, 1, "Not implemented yet");
  }  
    
  if(strcmp(pot_type, "single") == 0 &&
     strcmp(problem_type, "driv") == 0) {
    ierr = Driv1dCalc2_single_driv(energy); CHKERRQ(ierr);
  }

  if(strcmp(pot_type, "double") == 0 &&
     strcmp(problem_type, "driv") == 0) {
    ierr = Driv1dCalc2_double_driv(energy); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode Driv1dWriteWaveFunc() {
  PetscErrorCode ierr;
  if(use_func_view) {
    PrintTimeStamp(comm, "Writing", NULL);
    if(strcmp(pot_type, "single")==0) {
      ierr = FEMInfViewFunc(fem, c0, viewer); CHKERRQ(ierr);
    } else {
      ierr = FEMInfViewFunc(fem, c1, viewer); CHKERRQ(ierr);
    }
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
    printf("\n");
    printf("==== Calculation start ====\n");
    printf("energy : %f\n", energy);
    ierr = Driv1dCalc1(energy); CHKERRQ(ierr);
    ierr = Driv1dCalc2(energy); CHKERRQ(ierr);
  }
  ierr = Driv1dWriteWaveFunc(); CHKERRQ(ierr);
  ierr = Driv1dDestroy();       CHKERRQ(ierr);

  return 0;
}

