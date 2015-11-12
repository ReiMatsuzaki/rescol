#ifndef DVR_H
#define DVR_H

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <petscmat.h>
#include "mat.h"
#include "bps.h"

struct _p_DVR {
  MPI_Comm comm;
  int nq; // # of audrature in each element
  BPS bps; // breakpoints
  PetscReal *xs;
  PetscReal *ws;
  PetscInt num_basis;
  PetscReal *xs_basis; 
  Mat D2_R1LSMat;  // matrix of d^2/dx^2 for LS basis
  Mat R2_R1LSMat;  // matrix of r^{-2} for LS basis
  Mat T;   // Transformation matrix for LS basis to connected LS.
  Mat TT;  // Transpose of T
  Mat T2;  // Transformation matrix for LS R2 basis to connected LS R2.
  Mat T2T; // Transpose of T2;
};

typedef struct _p_DVR* DVR;

// ------- external functions ---------
PetscErrorCode MatCreateTransformMat(PetscReal *ws, int nq, int ne, 
				     MPI_Comm comm, Mat *B);
int NumDVR(int nq, int ne);

// ------- Basic Methods ---------
PetscErrorCode DVRCreate(MPI_Comm comm, DVR *p_self);
PetscErrorCode DVRDestroy(DVR *p_self);

PetscErrorCode DVRView(DVR self, PetscViewer v);

// ---- Accessor ----
PetscErrorCode DVRSetKnots(DVR self, int nq, BPS bps);
PetscErrorCode DVRSetUp(DVR self);
PetscErrorCode DVRSetFromOptions(DVR self);

PetscErrorCode DVRBasisPsi(DVR self, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode DVRGetSize(DVR self, int *n);

// ------- R1Mat/R2Mat ----------
PetscErrorCode DVRCreateR1Mat(DVR self, Mat *M);

PetscErrorCode DVRSR1Mat(DVR self, Mat M);
PetscErrorCode DVRD2R1Mat(DVR self, Mat M);
PetscErrorCode DVRR2invR1Mat(DVR self, Mat M);
PetscErrorCode DVRENR1Mat(DVR self, int q, PetscReal a, Mat M);
PetscErrorCode DVREER2Mat(DVR self, int q, Mat M);

// -------- LSR1Mat/LSR2Mat --------
PetscErrorCode DVRCreateR1LSMat(DVR self, Mat *M);
PetscErrorCode DVRCreateR2LSMat(DVR self, Mat *M);

PetscErrorCode DVRSR1LSMat(DVR self, Mat M);
PetscErrorCode DVRD2R1LSMat(DVR self, Mat M);
PetscErrorCode DVRR2invR1LSMat(DVR self, Mat M);
PetscErrorCode DVRENR1LSMat(DVR self, int q, PetscReal a, Mat M);
PetscErrorCode DVREER2LSMat(DVR self, int q, Mat M);

PetscErrorCode DVRLSMatToMat(DVR self, Mat A, MatReuse s, Mat *B);
PetscErrorCode DVRR2LSMatToR2Mat(DVR self, Mat A, Mat *B);

#ifdef __cplusplus
}
#endif
#endif
