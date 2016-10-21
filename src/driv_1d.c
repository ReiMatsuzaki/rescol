#include "../include/math.h"
#include "../include/range.h"
#include "../include/pot.h"
#include "../include/fem_inf.h"
#include "../include/viewerfunc.h"

static char help[] = "solve 1d driven equation.";

typedef struct {
  MPI_Comm  comm;
  PetscInt  L;    // angular quantum number
  PetscReal Z;    // charge
  PetscBool use_V;// True=> use V
  Pot       V;    // additional potential
  Pot       S;    // driven term

  FEMInf    fem;  // grid method
  Range     energy_range;  // calculation energy range
  PetscBool use_func_view; // True=>print function
  ViewerFunc viewer;

  KSP ksp;
  Vec cvec, svec;
  Mat lmat;
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

PetscErrorCode Driv1dPrintIn(Driv1d self) {
  PrintTimeStamp(self->comm, "PrintIn", NULL);
  return 0;
}
PetscErrorCode Driv1dCalc(Driv1d self) {
  PrintTimeStamp(self->comm, "Calc", NULL);
  return 0;
}
PetscErrorCode Driv1dPrintOut(Driv1d self) {
  PrintTimeStamp(self->comm, "PrintOut", NULL);
  return 0;
}
PetscErrorCode Driv1dDestroy(Driv1d self) {
  PrintTimeStamp(self->comm, "Destroy", NULL);
  return 0;
}

int main(int argc, char **args) {

  PetscErrorCode ierr;
  Driv1d driv1d;
  MPI_Comm comm = MPI_COMM_SELF;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  ierr = Driv1dCreate(comm, &driv1d); CHKERRQ(ierr);
  ierr = Driv1dSetFromOptions(driv1d); CHKERRQ(ierr);
  ierr = Driv1dPrintIn(driv1d);  CHKERRQ(ierr);
  ierr = Driv1dCalc(driv1d);  CHKERRQ(ierr);
  ierr = Driv1dPrintOut(driv1d);  CHKERRQ(ierr);
  ierr = Driv1dDestroy(driv1d);  CHKERRQ(ierr);
  return 0;
  
}


