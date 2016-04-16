#ifndef H_TEMPLATE_H
#define H_TEMPLATE_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <petscmat.h>
#include <rescol/oce1.h>

struct _p_H2plus {

  int nr;
  int ny;

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
typedef struct _p_H2plus* OceH2plus;
PetscErrorCode H2PlusMatMultH(Mat H, Vec x, Vec y);
PetscErrorCode H2PlusMatMultS(Mat S, Vec x, Vec y);

PetscErrorCode H2plusCalc(OceH2plus self, OCE1 oce);
PetscErrorCode H2PlusCreateMat(OceH2plus self, Mat *M);
PetscErrorCode H2plusHMat(OceH2plus self, Mat H);
PetscErrorCode H2plusSMat(OceH2plus self, Mat S, PetscBool *is_id);

#ifdef __cplusplus
}
#endif
#endif
