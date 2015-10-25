#ifndef FD_H
#define FD_H

#include <stdio.h>
#include <petscmat.h>
#include "mat.h"

struct _p_FD {
  MPI_Comm comm;
  PetscScalar h;
  PetscInt num;
};

typedef struct _p_FD* FD;

// ---- Basic Methods ----
PetscErrorCode FDCreate(FD *fd, int num_xs, double xmax, MPI_Comm comm);
PetscErrorCode FDCreateFromOptions(FD *fd, MPI_Comm comm);
PetscErrorCode FDDestroy(FD *fd);
PetscErrorCode FDFPrintf(FD this, FILE *file, int lvl);

// ---- Matrix ----
PetscErrorCode FDSetSR1Mat(FD this, Mat *M);
PetscErrorCode FDSetD2R1Mat(FD this, Mat *M);
PetscErrorCode FDSetR2invR1Mat(FD this, Mat *M);
PetscErrorCode FDSetENR1Mat(FD this, int q, PetscScalar a, Mat *M);
PetscErrorCode FDSetEER2Mat(FD this, int q, Mat *M);

// ---- Vector ----
PetscErrorCode FDGuessHEig(FD this, int n, int l, PetscScalar z, Vec *v);

#endif
