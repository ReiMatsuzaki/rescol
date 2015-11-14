#ifndef BSPLINE_H
#define BSPLINE_H

#ifdef __cplusplus
extern "C" {
#endif 
#include <stdio.h>
#include <petscmat.h>
#include "mat.h"
#include "bps.h"
#include "cscaling.h"

struct _p_BSS {
  MPI_Comm comm;
  int order;    // equal to number of quadrature in each ele
  BPS bps; // break points
  int num_ele;  // number of finite elements
  int num_basis; 
  int* b_idx_list; // basis index list;

  CScaling c_scaling; // object for complex scaling

  int num_ts; // number of ts_r and ts_s
  PetscReal* ts_r; // overlapped knots points
  PetscScalar* ts_s; // overlapped knots points
  
  PetscReal* xs; // quadrature points
  PetscScalar* xs_s; // quadrature points (Scalar)
  PetscReal* ws; // weight
  PetscScalar *qrs; // coordinate transform
  PetscScalar *Rrs; // integral of qrs

  PetscReal* vals; // bspline values on quadrature points
  PetscReal* derivs; // derivative values

  PetscBool set_knots; 
  PetscBool setup; 
};
typedef struct _p_BSS* BSS;

// ----- external functions -----
int NumBSpline(int order, int num_ele);
int HasNon0Value(int order, int i, int j);
PetscErrorCode CalcBSpline(int order, PetscReal* ts_r, PetscScalar* ts_s, 
			   int i, double x, PetscScalar* y);
PetscErrorCode CalcDerivBSpline(int order, double* ts, PetscScalar* ts_s, 
				int i, double x, PetscScalar* y);


PetscErrorCode Non0QuadIndex(int a, int c, int k, int nq, int* i0, int* i1);

// ---- Basic Methods ----
PetscErrorCode BSSCreate(MPI_Comm comm, BSS *p_self);
PetscErrorCode BSSDestroy(BSS *p_self);

PetscErrorCode BSSView(BSS self, PetscViewer v);
PetscErrorCode BSSCheck(BSS self);

PetscErrorCode BSSSetKnots(BSS self, int order, BPS bps);
PetscErrorCode BSSSetCScaling(BSS self, CScaling cscaling);
PetscErrorCode BSSSetUp(BSS self);
PetscErrorCode BSSSetFromOptions(BSS self);

PetscErrorCode BSSPsi(BSS self, Vec c, PetscReal x, PetscScalar *y);
PetscErrorCode BSSBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y);
PetscErrorCode BSSDerivBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y);
PetscErrorCode BSSGetSize(BSS self, int *n);
  

PetscErrorCode BSSCreateR1Mat(BSS self, Mat *M);
PetscErrorCode BSSCreateR2Mat(BSS self, Mat *M);
PetscErrorCode BSSCreateR1Vec(BSS self, Vec *v);

PetscErrorCode BSSSR1Mat(BSS self, Mat S);
PetscErrorCode BSSR2invR1Mat(BSS self, Mat M);
PetscErrorCode BSSD2R1Mat(BSS self, Mat D);
PetscErrorCode BSSENR1Mat(BSS self, int q, PetscReal  a, Mat V);
PetscErrorCode BSSPotR1Mat(BSS self, PF pot, Mat M);
PetscErrorCode BSSPotR1Vec(BSS self, PF pot, Vec v);
PetscErrorCode BSSEER2Mat(BSS self, int q, Mat V);
PetscErrorCode BSSEER2Mat_ver1(BSS self, int q, Mat V);

/*
PetscErrorCode BSSSetSR1Mat(BSS self, Mat *S);
PetscErrorCode BSSSetR2invR1Mat(BSS self, Mat *M);
PetscErrorCode BSSSetD2R1Mat(BSS self, Mat *D);
PetscErrorCode BSSSetENR1Mat(BSS self, int q, PetscScalar a, Mat *D);
PetscErrorCode BSSSetPotR1Mat(BSS self, POT pot, Mat *M);
PetscErrorCode BSSSetEER2Mat(BSS self, int q, Mat *V);
*/

#ifdef __cplusplus
}
#endif
#endif
