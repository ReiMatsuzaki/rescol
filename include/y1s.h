#ifndef Y1S_H
#define Y1S_H
#ifdef __cplusplus
extern "C" {
#endif 
#include "angmoment.h"

// ---- function ----
PetscReal Y1RedYq(int j1, int j2, int j3);
PetscReal Y1EleYqk(int j1, int j2, int j3, int m1, int m2, int m3);
PetscReal Y1ElePq(int j1, int q, int j2, int m1, int m2);

// ---- Set of Y1 ----
struct _p_Y1s {
  MPI_Comm comm;
  int num;
  int* ls;
  int m;  
};
typedef struct _p_Y1s* Y1s;
PetscErrorCode Y1sCreate(MPI_Comm comm, Y1s *p_self);
PetscErrorCode Y1sDestroy(Y1s *p_self);

PetscErrorCode Y1sView(Y1s y1s, PetscViewer viewer);

PetscErrorCode Y1sSet(Y1s self, int m, int g_or_u, int lmax);
PetscErrorCode Y1sSetOne(Y1s self, int m, int l);
PetscErrorCode Y1sSetFromOptions(Y1s self);

PetscErrorCode Y1sGetSize(Y1s y1s, int *n);
PetscErrorCode Y1sGetMaxL(Y1s ys, int *lmax);

PetscErrorCode Y1sCreateY1Mat(Y1s self, Mat *M);
PetscErrorCode Y1sCreateY1MatOther(Y1s self, Y1s other, Mat *M);
PetscErrorCode Y1sSY1Mat(Y1s self, Mat M);
PetscErrorCode Y1sLambdaY1Mat(Y1s self, Mat M);
PetscErrorCode Y1sPqY1Mat(Y1s self, int q, Mat M, PetscBool *non0);
PetscErrorCode Y1sYqkY1MatOther(Y1s self, Y1s other, int q, int k, Mat M, PetscBool *non0);

#ifdef __cplusplus
}
#endif 
#endif
