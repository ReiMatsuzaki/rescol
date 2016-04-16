#ifndef FD_H
#define FD_H
#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <petscmat.h>
#include "mat.h"

struct _p_FD {
  MPI_Comm comm;
  PetscReal h;
  PetscInt num;
};

typedef struct _p_FD* FD;

// ---- Basic Methods ----
PetscErrorCode FDCreate(MPI_Comm comm, FD *p_self);
PetscErrorCode FDDestroy(FD *p_self);

PetscErrorCode FDView(FD self, PetscViewer v);

// ---- Accessor ----
PetscErrorCode FDSetMesh(FD self, int num_xs, PetscReal xmax);
PetscErrorCode FDSetFromOptions(FD self);

PetscErrorCode FDGetSize(FD self, int *n);

// ---- Matrix/Vector ----
PetscErrorCode FDCreateR1Mat(FD self, Mat *M);
PetscErrorCode FDCreateR1Vec(FD self, Vec *M);
PetscErrorCode FDCreateR2Mat(FD self, Mat *M);

PetscErrorCode FDSR1Mat(FD self, Mat M);
PetscErrorCode FDD2R1Mat(FD self, Mat M);
PetscErrorCode FDR2invR1Mat(FD self, Mat M);
PetscErrorCode FDENR1Mat(FD self, int q, PetscScalar a, Mat M);
PetscErrorCode FDEER2Mat(FD self, int q, Mat M);

PetscErrorCode FDGuessHEig(FD self, int n, int l, PetscScalar z, Vec v);

#ifdef __cplusplus
}
#endif
#endif
