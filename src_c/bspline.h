#ifndef BSPLINE_H
#define BSPLINE_H

#include <stdio.h>
#include <petscmat.h>
#include "mat.h"
#include "bps.h"

struct _p_BSS {
  MPI_Comm comm;
  int order;    // equal to number of quadrature in each ele
  BPS bps;
  int num_ele;  // number of finite elements
  int num_basis; 
  PetscScalar rmax;
  int* b_idx_list; // basis index list;
  PetscScalar* ts; // overlapped knots points
  PetscScalar* xs; // quadrature points
  PetscScalar* ws; // weight
  PetscScalar* vals; // bspline values on quadrature points
  PetscScalar* derivs; // derivative values
};

typedef struct _p_BSS* BSS;

// ----- external functions -----
int NumBSpline(int order, int num_ele);
int HasNon0Value(int order, int i, int j);
PetscErrorCode CalcBSpline(int order, double* ts, int i, double x, double* y);
PetscErrorCode CalcDerivBSpline(int order, double* ts, int i, double x, double* y);

// ---- Basic Methods ----
PetscErrorCode BSSCreate(BSS *bss, int order, BPS bps, MPI_Comm comm);
PetscErrorCode BSSCreateFromOptions(BSS *bss, MPI_Comm comm);
PetscErrorCode BSSDestroy(BSS *bss);
PetscErrorCode BSSFPrintf(BSS this, FILE* file, int lvl);
PetscErrorCode BSSBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode BSSDerivBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y);

// ---- Matrix ----
PetscErrorCode BSSInitR1Mat(BSS this, Mat *M);
PetscErrorCode BSSInitR2Mat(BSS this, Mat *M);

PetscErrorCode BSSCalcSR1Mat(BSS this, Mat S, InsertMode mode);
PetscErrorCode BSSCalcR2invR1Mat(BSS this, Mat M, InsertMode mode);
PetscErrorCode BSSCalcD2R1Mat(BSS this, Mat D, InsertMode);  
PetscErrorCode BSSCalcENR1Mat(BSS this, int q, PetscScalar  a, Mat V, InsertMode); 
PetscErrorCode BSSCalcEER2Mat(BSS this, int q, Mat V, InsertMode);
PetscErrorCode BSSCalcEER2Mat_ver1(BSS this, int q, Mat V, InsertMode);

PetscErrorCode BSSSetSR1Mat(BSS this, Mat *S);
PetscErrorCode BSSSetR2invR1Mat(BSS this, Mat *M);
PetscErrorCode BSSSetD2R1Mat(BSS this, Mat *D);
PetscErrorCode BSSSetENR1Mat(BSS this, int q, PetscScalar a, Mat *D);
PetscErrorCode BSSSetEER2Mat(BSS this, int q, Mat *V);

#endif
