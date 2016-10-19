#ifndef POT_H
#define POT_H
#ifdef __cplusplus
extern "C" {
#endif
#include <petscpf.h>
typedef PF Pot;

int str_nele(const char *str, char key);
PetscErrorCode PotCreate(MPI_Comm comm, Pot *p_self);
PetscErrorCode PotSetFromOptions(Pot self);
PetscErrorCode PotSetFromStr(Pot self, const char str[]);
PetscErrorCode PotSetFromOptions2(Pot self, const char prefix[], PetscBool *find);
PetscErrorCode PotSetHarm(Pot self, PetscScalar a);
PetscErrorCode PotSetPower(Pot self, PetscScalar a, PetscInt n);
PetscErrorCode PotSetCoulombNE(Pot self, int q, PetscScalar a, PetscScalar zz);
PetscErrorCode PotSetSlater(Pot self, PetscScalar a, int n, PetscScalar z);
PetscErrorCode PotSetMorse(Pot self, PetscScalar D0, PetscScalar a, PetscScalar Re);
PetscErrorCode PotSetRBessel(Pot self, int L, double k);
PetscErrorCode PotSetCombination(Pot self, int num, Pot *pfs);
PetscErrorCode PotSetProduct(Pot self, int num, Pot *pfs);

#ifdef __cplusplus
}
#endif
#endif
