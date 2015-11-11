#ifndef POT_H
#define POT_H
#ifdef __cplusplus
extern "C" {
#endif
#include <petscpf.h>

PetscErrorCode PFSetPotentialFromOptions(PF self);
PetscErrorCode PFSetHarm(PF self, PetscScalar a);
PetscErrorCode PFSetPower(PF self, PetscScalar a, PetscInt n);
PetscErrorCode PFSetCoulomnNE(PF self, int q, PetscScalar a);
PetscErrorCode PFSetSlater(PF self, PetscScalar a, int n, PetscScalar z);
PetscErrorCode PFSetMorse(PF self, PetscScalar D0, PetscScalar a, PetscScalar Re);

#ifdef __cplusplus
}
#endif
#endif
