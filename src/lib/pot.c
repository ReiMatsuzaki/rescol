#include <rescol/pot.h>

// ---- Potential Object ----
PetscErrorCode POTCreate(POT *pot) {
  POT _pot;
  PetscNew(&_pot);
  *pot = _pot;
  return 0;
}
PetscErrorCode POTCalc(POT pot, PetscScalar x, PetscScalar *y) {
  *y = pot->Calc(x, pot->vs);
  return 0;
}
PetscErrorCode POTView(POT pot) {
  PetscErrorCode ierr;
  ierr = pot->View(pot->vs);
  return 0;
}
PetscErrorCode POTViewFunc(POT self, ViewerFunc viewer) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_WORLD;
  ierr = ViewerFuncCheckAcrive(viewer); CHKERRQ(ierr);

  PetscViewerType type;
  PetscViewerGetType(viewer->base, &type);
  if(strcmp(type, "ascii") != 0) {
    char msg[100]; sprintf(msg, "unsupported type: %s", type);
    SETERRQ(comm, 1, msg);
  }

  PetscInt num;
  PetscReal *xs; 
  ViewerFuncGetXs(viewer, &num, &xs);

  for(int i = 0; i < num; i++) {
    PetscScalar y;
    POTCalc(self, xs[i], &y);
    PetscReal re, im;
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(y);
    im  =PetscImaginaryPart(y);
#else
    re = y;
    im = 0.0;
#endif
    PetscViewerASCIIPrintf(viewer->base, "%lf %lf %lf\n", xs[i], re, im);
  }
  
  return 0;
}
PetscErrorCode POTDestroy(POT *pot) {
  PetscErrorCode ierr;
  
  ierr = PetscFree((*pot)->vs); CHKERRQ(ierr);
  ierr = PetscFree(*pot); CHKERRQ(ierr);
  return 0;
}
PetscBool POTIsType(POT pot, char *name) {
  return (strcmp(pot->name, name) == 0);
}

// ---- Harmonic Potential ----
PetscScalar POTHarmCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * 0.5 * x * x;
}
PetscErrorCode POTHarmView(PetscScalar *vs) {
  PetscPrintf(PETSC_COMM_SELF, ">>> POT harm >>>\n");  
  PetscPrintf(PETSC_COMM_SELF, "       a  2 \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "v(x) = - x  \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "       2    \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "a = %f\n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "<<< POT harm <<<\n");
  return 0;
}
PetscErrorCode POTHarmCreate(POT *pot, PetscScalar a) {
  
  POT self;
  POTCreate(&self);
  PetscScalar *vs; 
  PetscMalloc1(1, &vs); vs[0] = a;

  strcpy(self->name, "harm");
  self->num = 1;
  self->vs = vs;
  self->Calc = POTHarmCalc;
  self->View = POTHarmView;

  *pot = self;
  return 0;
}

// ---- Power Potential ----
PetscScalar POTPowerCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * pow(x, vs[1]);
}
PetscErrorCode POTPowerView(PetscScalar *vs) {
  PetscPrintf(PETSC_COMM_SELF, ">>> POT Power >>>\n");
  PetscPrintf(PETSC_COMM_SELF, "          n \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "v(x) = A x  \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "            \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "A = %f\n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "n = %f\n", vs[1]);
  PetscPrintf(PETSC_COMM_SELF, "<<< POT Power <<<\n");
  return 0;
}
PetscErrorCode POTPowerCreate(POT *pot, PetscScalar a, PetscScalar n) {

  POTCreate(pot);
  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = a;
  vs[1] = n;
  strcpy((*pot)->name, "power");
  (*pot)->num = 2;
  (*pot)->vs = vs;
  (*pot)->Calc = POTPowerCalc;
  (*pot)->View = POTPowerView;
  return 0;

}

// ---- Coulomb Potential ----
PetscScalar POTCoulombCalc(PetscScalar x, PetscScalar *vs) {

  PetscReal a_r = PetscRealPart(vs[1]);
  PetscReal x_r = PetscRealPart(x);
  PetscScalar g = a_r > x_r ? vs[1] : x;
  PetscScalar s = a_r < x_r ? vs[1] : x;
  return pow(s/g, vs[0])/g;
}
PetscErrorCode POTCoulombView(PetscScalar *vs) {
  PetscPrintf(PETSC_COMM_SELF, ">>> POT Coulomb >>>\n");
  PetscPrintf(PETSC_COMM_SELF, "          min[a,x]^q   \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "v(x) =  -------------- \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "        max[a,x]^{q+1} \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "q = %f\n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "a = %f\n", vs[1]);
  PetscPrintf(PETSC_COMM_SELF, "<<< POT Coulomb <<<\n");
  return 0;
}
PetscErrorCode POTCoulombCreate(POT *pot, PetscScalar q, PetscScalar a) {

  POTCreate(pot);
  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = q;
  vs[1] = a;
  strcpy((*pot)->name, "coulomb");
  (*pot)->num = 2;
  (*pot)->vs = vs;
  (*pot)->Calc = POTCoulombCalc;
  (*pot)->View = POTCoulombView;
  return 0;
}

// ---- Slater Potential ----
PetscScalar POTSlaterCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * x * x * exp(-vs[1]*x);
}
PetscErrorCode POTSlaterView(PetscScalar *vs) {
  PetscPrintf(PETSC_COMM_SELF, ">>> POT Slater >>>\n");
  PetscPrintf(PETSC_COMM_SELF, "             2  -zx \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "v(x) =  V_0 x  e    \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "V_0 = %f\n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "z   = %f\n", vs[1]);
  PetscPrintf(PETSC_COMM_SELF, "<<< POT Slater <<<\n");  
  return 0;
}
PetscErrorCode POTSlaterCreate(POT *pot, PetscScalar v0, PetscScalar z) {
  POTCreate(pot);
  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = v0;
  vs[1] = z;
  strcpy((*pot)->name, "slater");
  (*pot)->num = 2;
  (*pot)->vs = vs;
  (*pot)->Calc = POTSlaterCalc;
  (*pot)->View = POTSlaterView;
  return 0;
}

// ---- Morse Potential ----
PetscScalar POTMorseCalc(PetscScalar x, PetscScalar *vs) {
  PetscScalar D0 = vs[0];
  PetscScalar a = vs[1];
  PetscScalar Re = vs[2];
  PetscScalar t = 1.0 - exp(-a*(x-Re));
  return D0 * t * t;
}
PetscErrorCode POTMorseView(PetscScalar *vs) {
  PetscPrintf(PETSC_COMM_SELF, ">>> POT Morse >>>\n");
  PetscPrintf(PETSC_COMM_SELF, "                             2 \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "v(x) =  D_0 (1-exp[-a(x-Re)])  \n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "D_0 = %f\n", vs[0]);
  PetscPrintf(PETSC_COMM_SELF, "a   = %f\n", vs[1]);
  PetscPrintf(PETSC_COMM_SELF, "Re  = %f\n", vs[2]);
  PetscPrintf(PETSC_COMM_SELF, "<<< POT Morse <<<\n");  
  return 0;
}
PetscErrorCode POTMorseCreate(POT *p_self, PetscScalar D0, PetscScalar a, PetscScalar Re) {
  POTCreate(p_self);
  PetscScalar *vs; 
  PetscMalloc1(3, &vs); 
  vs[0] = D0;
  vs[1] = a;
  vs[2] = Re;
  strcpy((*p_self)->name, "morse");
  (*p_self)->num = 3;
  (*p_self)->vs = vs;
  (*p_self)->Calc = POTMorseCalc;
  (*p_self)->View = POTMorseView;
  return 0;
  return 0;
}

// ---- Create from options ----
PetscErrorCode POTCreateFromOptions(POT *pot, MPI_Comm comm) {
  
  char type[10];
  PetscErrorCode ierr;
  PetscBool find;

  ierr = PetscOptionsGetString(NULL, "-pot_type", type, 10, &find); CHKERRQ(ierr);

  if(!find) 
    SETERRQ(comm, 1, "options -pot_type is not found");
  
  if(strcmp(type, "harm") == 0) {
    PetscReal a;
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_a is not found");
    POTHarmCreate(pot, a);

  } else if(strcmp(type, "slater") == 0) {
    PetscReal v0, z;
    ierr = PetscOptionsGetReal(NULL, "-pot_v0", &v0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_v0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_z", &z, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_z is not found");
    POTSlaterCreate(pot, v0, z);

  } else if(strcmp(type, "morse") == 0) {
    PetscReal D0, a, Re;
    ierr = PetscOptionsGetReal(NULL, "-pot_D0", &D0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_D0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_z is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_Re", &Re, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_Re is not found");

    POTMorseCreate(pot, D0, a, Re);

  } else {
    char msg[100]; sprintf(msg, "unsupported pot_type: %s", type);
    SETERRQ(comm, 1, msg);
  }

  return 0;
}
