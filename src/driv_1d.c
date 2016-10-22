#include "../include/math.h"
#include "../include/range.h"
#include "../include/pot.h"
#include "../include/fem_inf.h"
#include "../include/viewerfunc.h"

static char help[] = "solve 1d driven equation:\n (E-T-V)psi=S\n";

typedef struct {
  // -- Inputs --
  MPI_Comm  comm;
  PetscInt  L;    // angular quantum number
  PetscReal Z;    // charge
  PetscBool use_V;// True=> use V
  Pot       V;    // additional potential
  Pot       S;    // driven term

  FEMInf     fem;  // grid method
  PetscBool  use_func_view; // True=>print function
  ViewerFunc viewer;
  Range      energy_range;  // calculation energy range

  KSP ksp;

  // -- Intermediate --
  //  Vec cvec, svec;
  //  Mat lmat;
  
  // -- Results --
  PetscScalar alpha;

} p_Driv1d;
typedef p_Driv1d* Driv1d;

PetscErrorCode Driv1dCreate(MPI_Comm comm, Driv1d *p_self) {

  PetscErrorCode ierr;
  
  PrintTimeStamp(comm, "Init", NULL);

  Driv1d self;
  ierr = PetscNew(&self); CHKERRQ(ierr);
  self->comm = comm;
  *p_self = self;

  self->L = 0;
  self->Z = 1.0;
  self->use_V = PETSC_FALSE;
  ierr = PotCreate(comm, &self->V); CHKERRQ(ierr);
  ierr = PotCreate(comm, &self->S); CHKERRQ(ierr);
  
  ierr = FEMInfCreate(comm, &self->fem); CHKERRQ(ierr);
  ierr = RangeCreate(comm, &self->energy_range); CHKERRQ(ierr);
  self->use_func_view = PETSC_FALSE;
  ierr = ViewerFuncCreate(comm, &self->viewer); CHKERRQ(ierr);
  
  ierr = KSPCreate(comm, &self->ksp); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode Driv1dSetFromOptions(Driv1d self) {

  PetscErrorCode ierr;
  PetscBool find;

  PrintTimeStamp(self->comm, "Set", NULL);

  PetscOptionsBegin(self->comm, "", "driv1d.c options", "none");

  // -- set L --
  ierr = PetscOptionsGetInt(NULL, NULL, "-L", &self->L, &find); CHKERRQ(ierr);
  if(!find)
    self->L = 0;
  if(self->L < 0)
    SETERRQ(self->comm, 1, "L must be zero or positive integer");

  // -- set Z --
  ierr = PetscOptionsGetReal(NULL, NULL, "-Z", &self->Z, &find); CHKERRQ(ierr);
  if(!find)
    self->Z = 0.0;


  // -- set V --
  ierr = PotSetFromOptions2(self->V, "V", &find); CHKERRQ(ierr);
  self->use_V = find;

  // -- set S --
  ierr = PotSetFromOptions2(self->S, "S", &find); CHKERRQ(ierr);
  if(!find)
    SETERRQ(self->comm, 1, "S is needed");

  // -- set FEM --
  ierr = FEMInfSetFromOptions(self->fem); CHKERRQ(ierr);

  // -- set energy range --
  ierr = RangeSetFromOptions(self->energy_range, "energy"); CHKERRQ(ierr);

  // -- set function viewer --
  ierr = ViewerFuncSetFromOptions(self->viewer, &find); CHKERRQ(ierr);
  self->use_func_view = find;

  PetscOptionsEnd();
  return 0;
}
PetscErrorCode Driv1dPrintIn(Driv1d self, PetscViewer viewer) {
  
  PrintTimeStamp(self->comm, "PrintIn", NULL);
  
  PetscViewerASCIIPrintf(viewer, "L:%d\n", self->L);
  
  PetscViewerASCIIPrintf(viewer, "L: %d\n", self->L);
  PetscViewerASCIIPrintf(viewer, "Z: %f\n", self->Z);
  PetscViewerASCIIPrintf(viewer, "use_V: %s\n", self->use_V ? "Yes" : "No");
  if(self->use_V) {
    PetscViewerASCIIPrintf(viewer, "V:\n");
    PetscViewerASCIIPushTab(viewer);
    PFView(self->V, viewer);
    PetscViewerASCIIPopTab(viewer);
  }
  PetscViewerASCIIPrintf(viewer, "S:\n");
  PetscViewerASCIIPushTab(viewer);
  PFView(self->S, viewer);
  PetscViewerASCIIPopTab(viewer);

  PetscViewerASCIIPrintf(viewer, "fem:\n");
  PetscViewerASCIIPushTab(viewer);
  FEMInfView(self->fem, viewer);
  PetscViewerASCIIPopTab(viewer);
  
  PetscViewerASCIIPrintf(viewer, "use_func_view: %s\n", self->use_func_view ? "Yes":"No");
  if(self->use_func_view) {
      PetscViewerASCIIPushTab(viewer);
      ViewerFuncView(self->viewer, viewer);
      PetscViewerASCIIPopTab(viewer);
  }
  
  PetscViewerASCIIPrintf(viewer, "energy_range: \n");
  PetscViewerASCIIPushTab(viewer);
  RangeView(self->energy_range, viewer);
  PetscViewerASCIIPopTab(viewer);
  
  return 0;
}
PetscErrorCode Driv1dCalc_mat(Driv1d self, Mat lmat, double energy) {

  PetscErrorCode ierr;

  int i[1] = {0}; 
  int j[1] = {0};
  PetscScalar v[1];
  
  ierr = FEMInfD2R1Mat(self->fem, lmat); CHKERRQ(ierr);
  ierr = MatScale(lmat, 0.5); CHKERRQ(ierr);
  MatGetValues(lmat, 1, i, 1, j, v);

  if(abs(self->Z) > 0.0000001) {
    Pot ZZ;
    Mat zmat;
    ierr = PotCreate(self->comm, &ZZ); CHKERRQ(ierr);
    ierr = FEMInfCreateMat(self->fem, 1, &zmat); CHKERRQ(ierr);

    ierr = PotSetPower(ZZ, -self->Z, -1); CHKERRQ(ierr);
    ierr = FEMInfPotR1Mat(self->fem, ZZ, zmat); CHKERRQ(ierr);
    ierr = MatAXPY(lmat, -1.0, zmat, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = PFDestroy(&ZZ); CHKERRQ(ierr);
    ierr = MatDestroy(&zmat); CHKERRQ(ierr);

    MatGetValues(lmat, 1, i, 1, j, v);
  }
  
  if(self->L != 0) {
    int L = self->L;
    Pot LL;
    Mat llmat;
    ierr = PotCreate(self->comm, &LL); CHKERRQ(ierr);
    ierr = FEMInfCreateMat(self->fem, 1, &llmat); CHKERRQ(ierr);
    
    ierr = PotSetPower(LL, L*(L+1)/2.0, -2); CHKERRQ(ierr);
    ierr = FEMInfPotR1Mat(self->fem, LL, llmat); CHKERRQ(ierr);
    ierr = MatAXPY(lmat, -1.0, llmat, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = PFDestroy(&LL); CHKERRQ(ierr);
    ierr = MatDestroy(&llmat); CHKERRQ(ierr);

    MatGetValues(lmat, 1, i, 1, j, v);
  }

  if(self->use_V) {
    Mat vmat;
    ierr = FEMInfCreateMat(self->fem, 1, &vmat); CHKERRQ(ierr);

    ierr = FEMInfPotR1Mat(self->fem, self->V, vmat); CHKERRQ(ierr);
    ierr = MatAXPY(lmat, -1.0, vmat, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

    ierr = MatDestroy(&vmat); CHKERRQ(ierr);
  }

  Mat smat;
  ierr = FEMInfCreateMat(self->fem, 1, &smat); CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(self->fem, smat); CHKERRQ(ierr);
  ierr = MatAXPY(lmat, +energy, smat, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

  MatGetValues(lmat, 1, i, 1, j, v);
  return 0;
}
PetscErrorCode Driv1dCalc(Driv1d self, PetscViewer viewer) {

  PetscErrorCode ierr;
  PrintTimeStamp(self->comm, "Calc", NULL);

  PetscReal energy;
  while(RangeNext(self->energy_range, &energy)) {
    PetscViewerASCIIPrintf(viewer, "energy = %f\n", energy);

    Mat lmat;
    ierr = FEMInfCreateMat(self->fem, 1, &lmat); CHKERRQ(ierr);
    ierr = Driv1dCalc_mat(self, lmat, energy);

    Vec svec;
    ierr = FEMInfCreateVec(self->fem, 1, &svec); CHKERRQ(ierr);
    ierr = FEMInfPotR1Vec(self->fem, self->S, svec); CHKERRQ(ierr);
  
    Vec cvec;
    ierr = FEMInfCreateVec(self->fem, 1, &cvec); CHKERRQ(ierr);
    ierr = KSPSetOperators(self->ksp, lmat, lmat); CHKERRQ(ierr);
    ierr = KSPSolve(self->ksp, svec, cvec); CHKERRQ(ierr);

    PetscScalar alpha;
    ierr = VecTDot(svec, cvec, &alpha); CHKERRQ(ierr);

    if(self->use_func_view) {
      ierr = FEMInfViewFunc(self->fem, cvec, self->viewer); CHKERRQ(ierr);
    }

    // -- destroy --
    ierr = MatDestroy(&lmat); CHKERRQ(ierr);
    ierr = VecDestroy(&svec); CHKERRQ(ierr);
    ierr = VecDestroy(&cvec); CHKERRQ(ierr);

    // -- print --
    PetscViewerASCIIPrintf(viewer, "alpha = (%f, %f)\n",
			   creal(alpha), cimag(alpha));
  }

  return 0;
}
PetscErrorCode Driv1dPrintOut(Driv1d self, PetscViewer viewer) {
  
  //  PrintTimeStamp(self->comm, "PrintOut", NULL);
  //  PetscViewerASCIIPrintf(viewer, "alpha = (%f, %f)\n",
  //			 creal(self->alpha),
  //			 cimag(self->alpha));
  //  return 0;
  
}
PetscErrorCode Driv1dDestroy(Driv1d self) {
  PrintTimeStamp(self->comm, "Destroy", NULL);
  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  Driv1d driv1d;
  MPI_Comm comm = MPI_COMM_SELF;

  // -- print compile time --
  printf("driv_1d program.\n");
  printf("Solve one dimensional driven equation:\n");
  printf("    (E-T-V)psi=S\n");
  printf("by grid method.");
  printf("Compile date: %s %s\n", __DATE__, __TIME__);
  
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  ierr = Driv1dCreate(comm, &driv1d); CHKERRQ(ierr);
  ierr = Driv1dSetFromOptions(driv1d); CHKERRQ(ierr);
  ierr = Driv1dPrintIn(driv1d, PETSC_VIEWER_STDOUT_SELF);  CHKERRQ(ierr);
  ierr = Driv1dCalc(driv1d, PETSC_VIEWER_STDOUT_SELF);  CHKERRQ(ierr);
  ierr = Driv1dPrintOut(driv1d, PETSC_VIEWER_STDOUT_SELF);  CHKERRQ(ierr);
  ierr = Driv1dDestroy(driv1d);  CHKERRQ(ierr);
  return 0;
  
}


