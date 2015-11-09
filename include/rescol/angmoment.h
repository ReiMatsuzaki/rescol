#ifndef ANGMOMENT_H
#define ANGMOMENT_H

#ifdef __cplusplus
extern "C" {
#endif 

#include <petscmat.h>
#include <gsl/gsl_sf_coupling.h>

#define SIGMA 0
#define PI 1
#define DELTA 2
#define PHI 3
#define GERADE 11
#define UNGERADE 12
#define PLUS 23
#define MINUS 24
#define ROT_SCALAR 101
#define ROT_VECTOR 102
typedef int RotSym;

// ---- utils ----
PetscReal wigner3j(int a, int b, int c, int d, int f, int e);
PetscReal wigner6j(int a, int b, int c, int d, int f, int e);
PetscBool TriangleQ(int j1, int j2, int j3);

// ---- function(y1) ----
PetscReal Y1RedYq(int j1, int j2, int j3);
PetscReal Y1EleYqk(int j1, int j2, int j3, int m1, int m2, int m3);
PetscReal Y1ElePq(int j1, int q, int j2, int m1, int m2);

// ---- Y1s ----
struct _p_Y1s {
  MPI_Comm comm;

  int num;
  int* ls;
  int m;  
};
typedef struct _p_Y1s* Y1s;
PetscErrorCode Y1sCreate(Y1s *y1s, MPI_Comm comm);
PetscErrorCode Y1sDestroy(Y1s *y1s);
PetscErrorCode Y1sSet(Y1s y1s, int m, int g_or_u, int lmax);
PetscErrorCode Y1sSetOne(Y1s y1s, int m, int l);
PetscErrorCode Y1sCreateFromOptions(Y1s *y1s, MPI_Comm comm);
PetscErrorCode Y1sView(Y1s y1s);
PetscErrorCode Y1sGetSize(Y1s y1s, int *n);
PetscErrorCode Y1sGetMaxL(Y1s ys, int *lmax);
PetscErrorCode Y1sCreateY1Mat(Y1s ys, Mat *M);
PetscErrorCode Y1sSetSY1Mat(Y1s self, Mat *M);
PetscErrorCode Y1sSetLambdaY1Mat(Y1s ys, Mat *M);
PetscErrorCode Y1sSetPqY1Mat(Y1s ys, int q, Mat *M);

// ---- function(y2) ----
typedef struct {
  int l1;
  int l2;
  int l;
  int m;
} Y2;
Y2 Y2Exchange(Y2 a);
PetscReal Y2EleYYq(Y2 a, int q, Y2 b);
PetscReal Y2ElePq12(Y2 a, int q, Y2 b);
PetscReal Y2ElePq1A(Y2 a, int q, Y2 b);
PetscReal Y2ElePq2A(Y2 a, int q, Y2 b);

// ---- Y2s ----
struct _p_Y2s {
  int num;
  Y2 *y2_list;
  MPI_Comm comm;
};
typedef struct _p_Y2s* Y2s;
PetscErrorCode Y2sCreate(Y2s *y2s, MPI_Comm comm);
PetscErrorCode Y2sCreateFromOptions(Y2s *y2s, MPI_Comm);
PetscErrorCode Y2sDestroy(Y2s *y2s);
PetscErrorCode Y2sSet(Y2s self, int m, int g_or_u, int p_or_m, int lmax);
PetscErrorCode Y2sView(Y2s self);
PetscErrorCode Y2sGetSize(Y2s self, int *n);
PetscErrorCode Y2sGetMaxL(Y2s self, int *L);
PetscErrorCode Y2sInitY2Mat(Y2s self, Mat *M);
PetscErrorCode Y2sSetSY2Mat(Y2s self, Mat *M);
PetscErrorCode Y2sSetGuessY2Vec(Y2s self, int L1, int L2, Vec *v);
PetscErrorCode Y2sSetLambda1Y2Mat(Y2s self, Mat *M);
PetscErrorCode Y2sSetLambda2Y2Mat(Y2s self, Mat *M);
PetscErrorCode Y2sSetPq1AY2Mat(Y2s self, int q, Mat *M);
PetscErrorCode Y2sSetPq2AY2Mat(Y2s self, int q, Mat *M);
PetscErrorCode Y2sSetPq12Y2Mat(Y2s self, int q, Mat *M);


#ifdef __cplusplus
}
#endif 
#endif
