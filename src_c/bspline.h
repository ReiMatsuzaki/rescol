#ifndef BSPLINE_H
#define BSPLINE_H

#include <stdio.h>
#include <petscmat.h>

struct _p_BSS {
  int order;    // equal to number of quadrature in each ele
  int num_ele;  // number of finite elements
  int num_basis; 
  int* b_idx_list; // basis index list;
  PetscScalar rmax;
  char knots_type[10];
  PetscScalar* ts; // overlapped knots points
  PetscScalar* zs; // non overlapped knots points
  PetscScalar* xs; // quadrature points
  PetscScalar* ws; // weight
  PetscScalar* vals; // bspline values on quadrature points
  PetscScalar* derivs; // derivative values
};

typedef struct _p_BSS* BSS;

int NumBSpline(int order, int num_ele);
int HasNon0Value(int order, int i, int j);
PetscErrorCode PartialCoulomb(int q, double r1, double r2, double *y);
PetscScalar LegGauss(int n, int i, PetscScalar* x, PetscScalar* w);
PetscErrorCode CalcBSpline(int order, double* ts, int i, double x, double* y);
PetscErrorCode CalcDerivBSpline(int order, double* ts, int i, double x, double* y);
PetscErrorCode CreateLinKnots(int num, double zmax, double *zs[]);

// Methods
PetscErrorCode BSSCreate(BSS *bss, int order, double* zs, int num_zs);
PetscErrorCode BSSCreateFromOptions(BSS *bss, MPI_Comm comm);
PetscErrorCode BSSDestroy(BSS *bss);
PetscErrorCode BSSFPrintf(BSS this, MPI_Comm comm, FILE* file, int lvl);
PetscErrorCode BSSBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode BSSDerivBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y);

PetscErrorCode BSSInitR1Mat(BSS this, MPI_Comm, Mat *M);
PetscErrorCode BSSInitR2Mat(BSS this, MPI_Comm, Mat *M);

PetscErrorCode BSSCalcSR1Mat(BSS this, Mat S, InsertMode mode);
PetscErrorCode BSSCalcR2invR1Mat(BSS this, Mat M, InsertMode mode);
PetscErrorCode BSSCalcD2R1Mat(BSS this, Mat D, InsertMode);  
PetscErrorCode BSSCalcENR1Mat(BSS this, int q, double a, Mat V, InsertMode); 
PetscErrorCode BSSCalcEER2Mat(BSS this, int q, Mat V, InsertMode);
PetscErrorCode BSSCalcEER2Mat_ver1(BSS this, int q, Mat V, InsertMode);

PetscErrorCode BSSSetSR1Mat(BSS this, MPI_Comm comm, Mat *S);
PetscErrorCode BSSSetR2invR1Mat(BSS this, MPI_Comm comm, Mat *M);
PetscErrorCode BSSSetD2R1Mat(BSS this, MPI_Comm comm, Mat *D);
PetscErrorCode BSSSetENR1Mat(BSS this, int q, double a, MPI_Comm comm, Mat *D);
PetscErrorCode BSSSetEER2Mat(BSS this, int q, MPI_Comm comm, Mat *V);
PetscErrorCode BSSSetUR1R2Mat(BSS this, Mat *U);
PetscErrorCode BSSSetEER2MatGreen(BSS this, int q, MPI_Comm, Mat *V);

#endif
