#ifndef OCE1_H
#define OCE1_H
#ifdef __cplusplus
extern "C" {
#endif 

/*
  matrix caluculator using one center expansion for physical system
*/

#include "fem_inf.h"
#include "y1s.h"
#include "synthesize.h"

struct _p_OCE1 {
  MPI_Comm comm;
  PetscReal mu;
  FEMInf fem;
  Y1s y1s;

  Mat s_r;
  Mat s_y;
};
typedef struct _p_OCE1* OCE1;

PetscErrorCode OCE1Create(MPI_Comm comm, OCE1 *p_self);
PetscErrorCode OCE1Destroy(OCE1 *p_self);

PetscErrorCode OCE1View(OCE1 self, PetscViewer v);
PetscErrorCode OCE1ViewFunc(OCE1 self, Vec c, ViewerFunc v);

PetscErrorCode OCE1Set(OCE1 self, FEMInf fem, Y1s y1s);
PetscErrorCode OCE1SetFromOptions(OCE1 self);

PetscErrorCode OCE1GetSizes(OCE1 self, int *n_r, int *n_y);

PetscErrorCode OCE1Fit(OCE1 self, PF pf, int L, KSP ksp, Vec c);
PetscErrorCode OCE1Psi(OCE1 self, Vec cs, int L, int M, Vec xs, Vec ys);
PetscErrorCode OCE1PsiOne(OCE1 self, Vec cs, int L, int M, PetscScalar x, PetscScalar *y);
PetscErrorCode OCE1CreateMat(OCE1 self, Mat *M);
PetscErrorCode OCE1CreateMatOther(OCE1 self, OCE1 other, Mat *M);
PetscErrorCode OCE1CreateVec(OCE1 self, Vec *v);
PetscErrorCode OCE1SMat(OCE1 self,MatReuse scall,  Mat *M, PetscBool *is_id);
  //PetscErrorCode OCE1SMatNullable(OCE1 self, Mat M);
PetscErrorCode OCE1TMat(OCE1 self, MatReuse scall, Mat *M);
PetscErrorCode OCE1PotMat(OCE1 self, RotSym sym, Pot pot, MatReuse scall, Mat *M);
PetscErrorCode OCE1PlusPotMat(OCE1 self, RotSym sym, Pot pot, Mat M);
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat M);
PetscErrorCode OCE1H2PlusMat(OCE1 self, PetscReal a, PetscReal z, Mat *H, Mat *S, PetscBool *is_id);
PetscErrorCode OCE1ZMat(OCE1 a, OCE1 b, MatReuse scall, Mat *M);

struct _p_OceH2plus {

  PetscReal a;
  PetscReal z;

  Mat s_r1;
  Mat d2_r1;
  Mat r2inv_r1;

  Mat s_y1;
  Mat lambda_y1;

  int nq;
  int *q;
  Mat *ne_r1;
  Mat *pq_y1;

};
typedef struct _p_OceH2plus* OceH2plus;
PetscErrorCode OCE1CreateH2plus(OCE1 self, PetscReal a, PetscReal z, OceH2plus *p_ctx);
PetscErrorCode OCE1H2plusMat(OCE1 self, OceH2plus ctx, Mat *H, Mat *S, PetscBool *is_id);
PetscErrorCode OCE1H2plusMat_direct(OCE1 self, OceH2plus ctx, Mat *H, Mat *S, PetscBool *is_id);

PetscErrorCode OCE1H2plusDestroy(OceH2plus *p_ctx);

  //PetscErrorCode OceH2plusMatMultH(Mat H, Vec x, Vec y);
  //PetscErrorCode OceH2plusMatMultS(Mat S, Vec x, Vec y);

/*
PetscErrorCode OceH2plusCreate(MPI_Comm comm, OceH2plus *p_self);
PetscErrorCode OceH2plusDestroy(OceH2plus *p_self);

PetscErrorCode OceH2plusCalc(OceH2plus self, FEMInf fem, Y1s y1s);
PetscErrorCode OceH2PlusCreateMat(OceH2plus self, Mat *M);
PetscErrorCode OceH2plusHMat(OceH2plus self, Mat H);
PetscErrorCode OceH2plusSMat(OceH2plus self, Mat S, PetscBool *is_id);
*/
#ifdef __cplusplus
}
#endif 
#endif
