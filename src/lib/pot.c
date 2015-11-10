#include <rescol/pot.h>

// ---- Potential Object ----
PetscErrorCode POTCreate(MPI_Comm comm, POT *p_self) {
  POT self;
  PetscNew(&self);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode POTDestroy(POT *p_self) {
  PetscErrorCode ierr;
  ierr = PetscFree((*p_self)->vs); CHKERRQ(ierr);
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POTView(POT self, PetscViewer v) {
  PetscErrorCode ierr;

  PetscViewerType type;
  PetscViewerGetType(v, &type);
  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "only ascii is supported");

  ierr = self->View(self->vs, v);
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

PetscErrorCode POTCalc(POT self, PetscScalar x, PetscScalar *y) {
  *y = self->Calc(x, self->vs);
  return 0;
}
PetscBool POTIsType(POT self, char *name) {
  return (strcmp(self->name, name) == 0);
}

// ---- Harmonic Potential ----
PetscScalar POTHarmCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * 0.5 * x * x;
}
PetscErrorCode POTHarmView(PetscScalar *vs, PetscViewer v ) {

  PetscViewerASCIIPrintf(v, ">>> POT harm >>>\n");  
  PetscViewerASCIIPrintf(v, "       a  2 \n", vs[0]);
  PetscViewerASCIIPrintf(v, "v(x) = - x  \n", vs[0]);
  PetscViewerASCIIPrintf(v, "       2    \n", vs[0]);
  PetscViewerASCIIPrintf(v, "a = %f\n", vs[0]);
  PetscViewerASCIIPrintf(v, "<<< POT harm <<<\n");
  return 0;
}
PetscErrorCode POTSetHarm(POT self, PetscScalar a) {

  PetscScalar *vs; 
  PetscMalloc1(1, &vs); vs[0] = a;

  strcpy(self->name, "harm");
  self->num = 1;
  self->vs = vs;
  self->Calc = POTHarmCalc;
  self->View = POTHarmView;

  return 0;
}

// ---- Power Potential ----
PetscScalar POTPowerCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * pow(x, vs[1]);
}
PetscErrorCode POTPowerView(PetscScalar *vs, PetscViewer v) {

  PetscViewerASCIIPrintf(v, ">>> POT Power >>>\n");
  PetscViewerASCIIPrintf(v, "          n \n", vs[0]);
  PetscViewerASCIIPrintf(v, "v(x) = A x  \n", vs[0]);
  PetscViewerASCIIPrintf(v, "            \n", vs[0]);
  PetscViewerASCIIPrintf(v, "A = %f\n", vs[0]);
  PetscViewerASCIIPrintf(v, "n = %f\n", vs[1]);
  PetscViewerASCIIPrintf(v, "<<< POT Power <<<\n");
  return 0;
}
PetscErrorCode POTSetPower(POT self, PetscScalar a, PetscScalar n) {

  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = a;
  vs[1] = n;
  strcpy(self->name, "power");
  self->num = 2;
  self->vs = vs;
  self->Calc = POTPowerCalc;
  self->View = POTPowerView;
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
PetscErrorCode POTCoulombView(PetscScalar *vs, PetscViewer v) {
  PetscViewerASCIIPrintf(v, ">>> POT Coulomb >>>\n");
  PetscViewerASCIIPrintf(v, "          min[a,x]^q   \n", vs[0]);
  PetscViewerASCIIPrintf(v, "v(x) =  -------------- \n", vs[0]);
  PetscViewerASCIIPrintf(v, "        max[a,x]^{q+1} \n", vs[0]);
  PetscViewerASCIIPrintf(v, "q = %f\n", vs[0]);
  PetscViewerASCIIPrintf(v, "a = %f\n", vs[1]);
  PetscViewerASCIIPrintf(v, "<<< POT Coulomb <<<\n");
  return 0;
}
PetscErrorCode POTSetCoulomb(POT self, PetscScalar q, PetscScalar a) {

  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = q;
  vs[1] = a;
  strcpy(self->name, "coulomb");
  self->num = 2;
  self->vs = vs;
  self->Calc = POTCoulombCalc;
  self->View = POTCoulombView;
  return 0;
}

// ---- Slater Potential ----
PetscScalar POTSlaterCalc(PetscScalar x, PetscScalar *vs) {
  return vs[0] * x * x * exp(-vs[1]*x);
}
PetscErrorCode POTSlaterView(PetscScalar *vs, PetscViewer v) {
  PetscViewerASCIIPrintf(v, ">>> POT Slater >>>\n");
  PetscViewerASCIIPrintf(v, "             2  -zx \n", vs[0]);
  PetscViewerASCIIPrintf(v, "v(x) =  V_0 x  e    \n", vs[0]);
  PetscViewerASCIIPrintf(v, "V_0 = %f\n", vs[0]);
  PetscViewerASCIIPrintf(v, "z   = %f\n", vs[1]);
  PetscViewerASCIIPrintf(v, "<<< POT Slater <<<\n");  
  return 0;
}
PetscErrorCode POTSetSlater(POT self, PetscScalar v0, PetscScalar z) {

  PetscScalar *vs; 
  PetscMalloc1(2, &vs); 
  vs[0] = v0;
  vs[1] = z;
  strcpy(self->name, "slater");
  self->num = 2;
  self->vs = vs;
  self->Calc = POTSlaterCalc;
  self->View = POTSlaterView;
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
PetscErrorCode POTMorseView(PetscScalar *vs, PetscViewer v) {
  PetscViewerASCIIPrintf(v, ">>> POT Morse >>>\n");
  PetscViewerASCIIPrintf(v, "                             2 \n", vs[0]);
  PetscViewerASCIIPrintf(v, "v(x) =  D_0 (1-exp[-a(x-Re)])  \n", vs[0]);
  PetscViewerASCIIPrintf(v, "D_0 = %f\n", vs[0]);
  PetscViewerASCIIPrintf(v, "a   = %f\n", vs[1]);
  PetscViewerASCIIPrintf(v, "Re  = %f\n", vs[2]);
  PetscViewerASCIIPrintf(v, "<<< POT Morse <<<\n");  
  return 0;
}
PetscErrorCode POTSetMorse(POT self, PetscScalar D0, PetscScalar a, PetscScalar Re) {

  PetscScalar *vs; 
  PetscMalloc1(3, &vs); 
  vs[0] = D0;
  vs[1] = a;
  vs[2] = Re;
  strcpy(self->name, "morse");
  self->num = 3;
  self->vs = vs;
  self->Calc = POTMorseCalc;
  self->View = POTMorseView;
  return 0;
}

// ---- Create from options ----
PetscErrorCode POTSetFromOptions(POT self) {
  
  char type[10];
  PetscErrorCode ierr;
  PetscBool find;

  ierr = PetscOptionsGetString(NULL, "-pot_type", type, 10, &find); CHKERRQ(ierr);

  if(!find) 
    SETERRQ(self->comm, 1, "options -pot_type is not found");
  
  if(strcmp(type, "harm") == 0) {
    PetscReal a;
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_a is not found");
    POTSetHarm(self, a);

  } else if(strcmp(type, "slater") == 0) {
    PetscReal v0, z;
    ierr = PetscOptionsGetReal(NULL, "-pot_v0", &v0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_v0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_z", &z, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_z is not found");
    POTSetSlater(self, v0, z);

  } else if(strcmp(type, "morse") == 0) {
    PetscReal D0, a, Re;
    ierr = PetscOptionsGetReal(NULL, "-pot_D0", &D0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_D0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_z is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_Re", &Re, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(self->comm, 1, "-pot_Re is not found");

    POTSetMorse(self, D0, a, Re);

  } else {
    char msg[100]; sprintf(msg, "unsupported pot_type: %s", type);
    SETERRQ(self->comm, 1, msg);
  }

  return 0;
}
