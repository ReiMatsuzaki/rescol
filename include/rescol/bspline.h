#ifndef BSPLINE_H
#define BSPLINE_H

#ifdef __cplusplus
extern "C" {
#endif 
#include <stdio.h>
#include <petscmat.h>
#include "mat.h"
#include "bps.h"
#include "pot.h"
#include "scale.h"

struct _p_BSS {
  MPI_Comm comm;
  int order;    // equal to number of quadrature in each ele
  BPS bps; // break points
  Scaler scaler; // object for complex scaling
  int num_ele;  // number of finite elements
  int num_basis; 
  //  PetscReal rmax;
  int* b_idx_list; // basis index list;
  
  PetscReal* ts_real; // overlapped knots points in real axix
  PetscScalar* ts; // overlapped knots points 
  PetscReal* xs; // quadrature points
  PetscReal* ws; // weight
  PetscScalar* vals; // bspline values on quadrature points
  PetscScalar* derivs; // derivative values
  PetscScalar *qrs; // coordinate transform
  PetscScalar *Rrs; // integral of qrs
};
typedef struct _p_BSS* BSS;

// ----- external functions -----
int NumBSpline(int order, int num_ele);
int HasNon0Value(int order, int i, int j);
PetscErrorCode CalcBSpline(int order, PetscReal* ts, int i, 
			   PetscReal x, PetscReal* y);
PetscErrorCode CalcDerivBSpline(int order, PetscReal* ts, int i, 
				PetscReal x, PetscReal* y);
PetscErrorCode Non0QuadIndex(int a, int c, int k, int nq, int* i0, int* i1);

// ---- Basic Methods ----
PetscErrorCode BSSCreate(BSS *bss, int order, BPS bps, Scaler scaler, 
			 MPI_Comm comm);
PetscErrorCode BSSCreateFromOptions(BSS *bss, MPI_Comm comm);
PetscErrorCode BSSDestroy(BSS *bss);
PetscErrorCode BSSFPrintf(BSS self, FILE* file, int lvl);
PetscErrorCode BSSBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y);
PetscErrorCode BSSDerivBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y);

// ---- Accessor ----
PetscErrorCode BSSGetSize(BSS self, int *n);

// ---- Matrix ----
PetscErrorCode BSSInitR1Mat(BSS self, Mat *M);
PetscErrorCode BSSInitR2Mat(BSS self, Mat *M);

PetscErrorCode BSSCalcSR1Mat(BSS self, Mat S, InsertMode mode);
PetscErrorCode BSSCalcR2invR1Mat(BSS self, Mat M, InsertMode mode);
PetscErrorCode BSSCalcD2R1Mat(BSS self, Mat D, InsertMode);  
PetscErrorCode BSSCalcENR1Mat(BSS self, int q, PetscReal  a, Mat V, InsertMode); 
PetscErrorCode BSSCalcEER2Mat(BSS self, int q, Mat V, InsertMode);
PetscErrorCode BSSCalcEER2Mat_ver1(BSS self, int q, Mat V, InsertMode);

PetscErrorCode BSSSetSR1Mat(BSS self, Mat *S);
PetscErrorCode BSSSetR2invR1Mat(BSS self, Mat *M);
PetscErrorCode BSSSetD2R1Mat(BSS self, Mat *D);
PetscErrorCode BSSSetENR1Mat(BSS self, int q, PetscScalar a, Mat *D);
PetscErrorCode BSSSetPotR1Mat(BSS self, POT pot, Mat *M);
PetscErrorCode BSSSetEER2Mat(BSS self, int q, Mat *V);

#ifdef __cplusplus
}
#endif
#endif
