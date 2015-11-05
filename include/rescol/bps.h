#ifndef BPS_H
#define BPS_H
#ifdef __cplusplus
extern "C" {
#endif

#include <petscmat.h>

/*
  Break points for finite elements methods
*/

struct _p_BPS {
  char type[10];
  PetscInt num_zs;
  PetscScalar *zs;
  MPI_Comm comm;
};

typedef struct _p_BPS* BPS;

// ------- Basic ----------
PetscErrorCode BPSCreate(BPS *bps, MPI_Comm comm);
PetscErrorCode BPSSetExp(BPS self, PetscScalar zmax, PetscInt num_zs, PetscScalar gamma);
PetscErrorCode BPSSetLine(BPS self, PetscScalar zmax, PetscInt num_zs);
PetscErrorCode BPSSetFromOptions(BPS self);
PetscErrorCode BPSDestroy(BPS *bpd);
PetscErrorCode BPSCheckPreallocated(BPS self);
PetscErrorCode BPSFPrintf(BPS self, FILE *file, int lvl);

// ------ Getter ---------
PetscErrorCode BPSGetZs(BPS self, PetscScalar **zs, PetscInt *num_zs);
PetscErrorCode BPSGetNumEle(BPS self, PetscInt *num_ele);
PetscErrorCode BPSGetZMax(BPS self, PetscScalar *zmax);

#ifdef __cplusplus
}
#endif
#endif