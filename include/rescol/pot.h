#ifndef POT_H
#define POT_H
#ifdef __cplusplus
extern "C" {
#endif
#include <petscpf.h>
typedef PF Pot;

PetscErrorCode PotCreate(MPI_Comm comm, Pot *p_self);
PetscErrorCode PotSetFromOptions(Pot self);
PetscErrorCode PotSetHarm(Pot self, PetscScalar a);
PetscErrorCode PotSetPower(Pot self, PetscScalar a, PetscInt n);
PetscErrorCode PotSetCoulombNE(Pot self, int q, PetscScalar a, PetscScalar zz);
PetscErrorCode PotSetSlater(Pot self, PetscScalar a, int n, PetscScalar z);
PetscErrorCode PotSetMorse(Pot self, PetscScalar D0, PetscScalar a, PetscScalar Re);

#ifdef __cplusplus
}
#endif
#endif
