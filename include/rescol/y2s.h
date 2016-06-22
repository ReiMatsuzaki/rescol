#ifndef Y2S_H
#define Y2S_H
#ifdef __cplusplus
extern "C" {
#endif 
#include "angmoment.h"

// ---- Coupled Y----
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

// ---- Set of Y2 ----
struct _p_Y2s {
  MPI_Comm comm;
  int num;
  Y2 *y2_list;  
};
typedef struct _p_Y2s* Y2s;
PetscErrorCode Y2sCreate(MPI_Comm comm, Y2s *p_self);
PetscErrorCode Y2sDestroy(Y2s *p_self);

PetscErrorCode Y2sView(Y2s self, PetscViewer viewer);

PetscErrorCode Y2sSetLM(Y2s self, int L, int M, int lmax);
PetscErrorCode Y2sSet(Y2s self, int m, int g_or_u, int p_or_m, int lmax);
PetscErrorCode Y2sSetFromOptions(Y2s p_self);

PetscErrorCode Y2sGetSize(Y2s self, int *n);
PetscErrorCode Y2sGetMaxL(Y2s self, int *L);

PetscErrorCode Y2sCreateY2Mat(Y2s self, Mat *M);
PetscErrorCode Y2sCreateY2Vec(Y2s self, Vec *V);
PetscErrorCode Y2sSY2Mat(Y2s self, Mat M);
PetscErrorCode Y2sGuessY2Vec(Y2s self, int L1, int L2, Vec v);
PetscErrorCode Y2sLambda1Y2Mat(Y2s self, Mat M, PetscBool *non0);
PetscErrorCode Y2sLambda2Y2Mat(Y2s self, Mat M, PetscBool *non0);
PetscErrorCode Y2sPq1AY2Mat(Y2s self, int q, Mat M, PetscBool *non0);
PetscErrorCode Y2sPq2AY2Mat(Y2s self, int q, Mat M, PetscBool *non0);
PetscErrorCode Y2sPq12Y2Mat(Y2s self, int q, Mat M, PetscBool *non0);

#ifdef __cplusplus
}
#endif
#endif
