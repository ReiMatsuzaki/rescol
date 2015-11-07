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
PetscErrorCode POTDestroy(POT *pot) {
  PetscErrorCode ierr;
  
  ierr = PetscFree((*pot)->vs); CHKERRQ(ierr);
  ierr = PetscFree(*pot); CHKERRQ(ierr);
  return 0;
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
  
  POTCreate(pot);
  PetscScalar *vs; 
  PetscMalloc1(1, &vs); vs[0] = a;
  (*pot)->num = 1;
  (*pot)->vs = vs;
  (*pot)->Calc = POTHarmCalc;
  (*pot)->View = POTHarmView;

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
  (*pot)->num = 2;
  (*pot)->vs = vs;
  (*pot)->Calc = POTSlaterCalc;
  (*pot)->View = POTSlaterView;
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

  } else if(strcmp(type, "slater")) {
    PetscReal v0, z;
    ierr = PetscOptionsGetReal(NULL, "-pot_v0", &v0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_v0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_z", &z, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_z is not found");
    POTSlaterCreate(pot, v0, z);

  } else {
    char msg[100]; sprintf(msg, "unsupported pot_type: %s", type);
    SETERRQ(comm, 1, msg);
  }

  return 0;
}
