#ifndef POT_H
#define POT_H

#ifdef __cplusplus
extern "C" {
#endif 
#include <petscmat.h>
#include <rescol/viewerfunc.h>

/*
  describing potential 
*/

// ---- Potential Class ----
struct _p_POT {
  MPI_Comm comm;
  char name[20];
  int num;
  PetscScalar *vs;
  PetscScalar (*Calc)(PetscScalar, PetscScalar*);  
  PetscErrorCode (*View)(PetscScalar*, PetscViewer);
};
typedef struct _p_POT* POT;

PetscErrorCode POTCreate(MPI_Comm comm, POT *p_self);
PetscErrorCode POTDestroy(POT *p_self);

PetscErrorCode POTView(POT self, PetscViewer v);
PetscErrorCode POTViewFunc(POT self, ViewerFunc viewer);

PetscErrorCode POTCalc(POT self, PetscScalar x, PetscScalar *y);
PetscBool POTIsType(POT pot, char *name);

// ---- Harmonic Potential ----
PetscErrorCode POTSetHarm(POT self, PetscScalar a);

// ---- Power Potential ----
PetscErrorCode POTSetPower(POT self, PetscScalar a, PetscScalar n);

// ---- Coulomb Potential ----
PetscErrorCode POTSetCoulomb(POT self, PetscScalar q, PetscScalar a);

// ---- Slater Potential ----
PetscErrorCode POTSetSlater(POT self, PetscScalar v0, PetscScalar z);

// ---- Mourse ----
PetscErrorCode POTSetMorse(POT self, PetscScalar D0, PetscScalar a, PetscScalar Re);

// ---- Create from options ----
PetscErrorCode POTSetFromOptions(POT self);

#ifdef __cplusplus
}
#endif

#endif
