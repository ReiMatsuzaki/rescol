#ifndef DVR_H
#define DVR_H
#ifdef __cplusplus
extern "C" {
#endif 

#include <stdio.h>
#include <petscmat.h>
#include "mat.h"
#include "bps.h"
#include "pot.h"
  //#include "cscaling.h"

struct _p_DVR {
  MPI_Comm comm;
  int nq; // # of audrature in each element
  BPS bps; // breakpoints
  PetscBool use_cscaling; // if true, use complex scaling 
  double R0;         // used in complex scaling
  double theta;      // used in scaling angle
  //  CScaling c_scaling; // object for complex scaling

  PetscInt num_basis;
  PetscReal *xs;
  PetscReal *xs_basis; 
  PetscScalar *xs_c; // used in complex scaling
  PetscScalar *ws_c; // used in complex scaling
  PetscScalar *xs_basis_c;  // used in complex scaling
  PetscScalar *ws_basis_c; // weight 
  
  Mat D2_R1LSMat;  // matrix of d^2/dx^2 for LS basis
  Mat R2_R1LSMat;  // matrix of r^{-2} for LS basis
  Mat T;   // Transformation matrix for LS basis to connected LS.
  Mat TT;  // Transpose of T
  Mat T2;  // Transformation matrix for LS R2 basis to connected LS R2.
  Mat T2T; // Transpose of T2;
};

typedef struct _p_DVR* DVR;

// ---- external functions ----
PetscErrorCode MatCreateTransformMat(PetscScalar *ws_c, int nq, int ne, 
				     MPI_Comm comm, Mat *B);
int NumDVR(int nq, int ne);
PetscErrorCode ValueLS(BPS bps, PetscScalar *xs_c,
		       PetscInt nq, 
		       PetscInt i, PetscInt m, PetscReal x,
		       PetscScalar x_c, PetscScalar *v, PetscBool *zeroq);

PetscErrorCode DerivLS(PetscScalar *xs_c, PetscScalar *ws_c, 
		       PetscInt ne, PetscInt nq, 
		       PetscInt i, PetscInt m, PetscInt mp, PetscScalar *v);

// ---- Basic Methods ----
PetscErrorCode DVRCreate(MPI_Comm comm, DVR *p_self);
PetscErrorCode DVRDestroy(DVR *p_self);

PetscErrorCode DVRView(DVR self, PetscViewer v);

// ---- Accessor ----
PetscErrorCode DVRSetKnots(DVR self, int nq, BPS bps);
PetscErrorCode DVRSetCScaling(DVR self, double R0, double theta);
PetscErrorCode DVRSetUp(DVR self);
PetscErrorCode DVRSetFromOptions(DVR self);

PetscErrorCode DVRBasisPsi(DVR self, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode DVRGetLSSize(DVR self, int *n);
PetscErrorCode DVRGetSize(DVR self, int *n);

// ----- Calculation ----
PetscErrorCode DVRPsiOne(DVR self, Vec cs, PetscReal x, PetscScalar *y);
PetscErrorCode DVRPsi(DVR self, Vec cs, Vec x, Vec y);
PetscErrorCode DVRDerivPsiOne(DVR self, Vec c, PetscReal x, PetscScalar *y);
PetscErrorCode DVRDerivPsi(DVR self, Vec c, Vec x, Vec y);


// ------- R1Mat/R2Mat ----------
PetscErrorCode DVRCreateR1Vec(DVR self, Vec *m);
PetscErrorCode DVRCreateR1Mat(DVR self, Mat *M);

PetscErrorCode DVRPotR1Vec(DVR self, Pot pot, Vec v);
PetscErrorCode DVRSR1Mat(DVR self, Mat M);
PetscErrorCode DVRD2R1Mat(DVR self, Mat M);
PetscErrorCode DVRPotR1Mat(DVR self, Pot pot, Mat M);
PetscErrorCode DVRR2invR1Mat(DVR self, Mat M);
PetscErrorCode DVRENR1Mat(DVR self, int q, PetscReal a, Mat M);
PetscErrorCode DVREER2Mat(DVR self, int q, Mat M);

// -------- LSR1Mat/LSR2Mat --------
PetscErrorCode DVRCreateR1LSVec(DVR self, Vec *m);
PetscErrorCode DVRCreateR1LSMat(DVR self, Mat *M);
PetscErrorCode DVRCreateR2LSMat(DVR self, Mat *M);

PetscErrorCode DVRPotR1LSVec(DVR self, Pot pot, Vec v);
PetscErrorCode DVRSR1LSMat(DVR self, Mat M);
PetscErrorCode DVRD2R1LSMat(DVR self, Mat M);
PetscErrorCode DVRR2invR1LSMat(DVR self, Mat M);
PetscErrorCode DVRENR1LSMat(DVR self, int q, PetscReal a, Mat M);
PetscErrorCode DVREER2LSMat(DVR self, int q, Mat M);

PetscErrorCode DVRLSVecToVec(DVR self, Vec A, Vec B);
PetscErrorCode DVRLSMatToMat(DVR self, Mat A, MatReuse s, Mat *B);
PetscErrorCode DVRR2LSMatToR2Mat(DVR self, Mat A, Mat *B);

#ifdef __cplusplus
}
#endif
#endif
