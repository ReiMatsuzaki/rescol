#include "pot.h"

// ---- Potential Object ----
PetscErrorCode POTCreate(POT *pot) {
  POT _pot;
  PetscNew(&_pot);
  *pot = _pot;
  return 0;
}
PetscScalar POTCalc(POT pot, PetscScalar x) {
  return pot->Calc(x, pot->vs);
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
  PetscScalar g = vs[1] > x ? vs[1] : x;
  PetscScalar s = vs[1] < x ? vs[1] : x;
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
