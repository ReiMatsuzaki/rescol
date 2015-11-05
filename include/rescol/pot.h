#ifndef POT_H
#define POT_H

#ifdef __cplusplus
extern "C" {
#endif 
#include <petscmat.h>

/*
  describing potential 
*/

// ---- Potential Class ----
struct _p_POT {
  int num;
  PetscScalar *vs;
  PetscScalar (*Calc)(PetscScalar, PetscScalar*);
  PetscErrorCode (*View)(PetscScalar*);
};
typedef struct _p_POT* POT;
PetscScalar POTCalc(POT pot, PetscScalar x);
PetscErrorCode POTView(POT pot);
PetscErrorCode POTDestroy(POT *pot);

// ---- Harmonic Potential ----
PetscErrorCode POTHarmCreate(POT *pot, PetscScalar a);

// ---- Power Potential ----
PetscErrorCode POTPowerCreate(POT *pot, PetscScalar a, PetscScalar n);

// ---- Coulomb Potential ----
PetscErrorCode POTCoulombCreate(POT *pot, PetscScalar q, PetscScalar a);

// ---- Slater Potential ----
PetscErrorCode POTSlaterCreate(POT *pot, PetscScalar v0, PetscScalar z);

#ifdef __cplusplus
}
#endif

#endif
