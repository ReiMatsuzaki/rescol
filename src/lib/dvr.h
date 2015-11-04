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
  PetscScalar *xs;
  PetscScalar *ws;
  PetscInt num_basis;
  PetscScalar *xs_basis; 
  Mat D2_R1LSMat;  // matrix of d^2/dx^2 for LS basis
  Mat R2_R1LSMat;  // matrix of r^{-2} for LS basis
  Mat T;   // Transformation matrix for LS basis to connected LS.
  Mat TT;  // Transpose of T
  Mat T2;  // Transformation matrix for LS R2 basis to connected LS R2.
  Mat T2T; // Transpose of T2;
};

typedef struct _p_DVR* DVR;

// ------- external functions ---------
PetscErrorCode MatCreateTransformMat(PetscScalar *ws, int nq, int ne, 
				     MPI_Comm comm, Mat *B);
int NumDVR(int nq, int ne);

// ------- Basic Methods ---------
PetscErrorCode DVRCreate(DVR *dvr, int nq, BPS bps, MPI_Comm comm);
PetscErrorCode DVRCreateFromOptions(DVR *dvr, MPI_Comm comm);
PetscErrorCode DVRDestroy(DVR *dvr);
PetscErrorCode DVRFPrintf(DVR self, FILE *file, int lvl);
PetscErrorCode DVRBasisPsi(DVR self, int i, PetscScalar x, PetscScalar *y);

// ---- Accessor ----
PetscErrorCode DVRGetSize(DVR self, int *n);

// ------- R1Mat ----------
PetscErrorCode DVRInitR1Mat(DVR self, Mat *M);
PetscErrorCode DVRSetSR1Mat(DVR self, Mat *M);
PetscErrorCode DVRSetD2R1Mat(DVR self, Mat *M);
PetscErrorCode DVRSetR2invR1Mat(DVR self, Mat *M);
PetscErrorCode DVRSetENR1Mat(DVR self, int q, double a, Mat *M);

// ------ R2Mat ---------
PetscErrorCode DVRSetEER2Mat(DVR self, int q, Mat *M);

// -------- LSR1Mat --------
PetscErrorCode DVRInitR1LSMat(DVR self, Mat *M);
PetscErrorCode DVRSetSR1LSMat(DVR self, Mat *M);
PetscErrorCode DVRSetD2R1LSMat(DVR self, Mat *M);
PetscErrorCode DVRSetR2invR1LSMat(DVR self, Mat *M);
PetscErrorCode DVRSetENR1LSMat(DVR self, int q, double a, Mat *M);

PetscErrorCode DVRLSMatToMat(DVR self, Mat A, Mat *B);

// ------- R2LSMat -----------
PetscErrorCode DVRInitR2LSMat(DVR self, Mat *M);
PetscErrorCode DVRSetEER2LSMat(DVR self, int q, Mat *M);

PetscErrorCode DVRR2LSMatToR2Mat(DVR self, Mat A, Mat *B);

#ifdef __cplusplus
}
#endif
#endif
