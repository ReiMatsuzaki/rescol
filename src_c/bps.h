#ifndef BPS_H
#define BPS_H

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
PetscErrorCode BPSSetExp(BPS this, PetscScalar zmax, PetscInt num_zs, PetscScalar gamma);
PetscErrorCode BPSSetLine(BPS this, PetscScalar zmax, PetscInt num_zs);
PetscErrorCode BPSSetFromOptions(BPS this);
PetscErrorCode BPSDestroy(BPS *bpd);
PetscErrorCode BPSCheckPreallocated(BPS this);
PetscErrorCode BPSFPrintf(BPS this, FILE *file, int lvl);

// ------ Getter ---------
PetscErrorCode BPSGetZs(BPS this, PetscScalar **zs, PetscInt *num_zs);
PetscErrorCode BPSGetNumEle(BPS this, PetscInt *num_ele);
PetscErrorCode BPSGetZMax(BPS this, PetscScalar *zmax);

#endif
