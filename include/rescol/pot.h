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
  char name[20];
  int num;
  PetscScalar *vs;
  PetscScalar (*Calc)(PetscScalar, PetscScalar*);  
  PetscErrorCode (*View)(PetscScalar*);
};
typedef struct _p_POT* POT;
PetscErrorCode POTCalc(POT pot, PetscScalar x, PetscScalar *y);
PetscErrorCode POTView(POT pot);
PetscErrorCode POTDestroy(POT *pot);
PetscBool POTIsType(POT pot, char *name);

// ---- Harmonic Potential ----
PetscErrorCode POTHarmCreate(POT *pot, PetscScalar a);

// ---- Power Potential ----
PetscErrorCode POTPowerCreate(POT *pot, PetscScalar a, PetscScalar n);

// ---- Coulomb Potential ----
PetscErrorCode POTCoulombCreate(POT *pot, PetscScalar q, PetscScalar a);

// ---- Slater Potential ----
PetscErrorCode POTSlaterCreate(POT *pot, PetscScalar v0, PetscScalar z);

// ---- Mourse ----
PetscErrorCode POTMorse(POT *p_self, PetscScalar D0, PetscScalar a, PetscScalar Re);

// ---- Create from options ----
PetscErrorCode POTCreateFromOptions(POT *pot, MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif
