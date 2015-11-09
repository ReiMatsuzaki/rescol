#ifndef OCE1_H
#define OCE1_H
#ifdef __cplusplus
extern "C" {
#endif 

/*
  matrix caluculator using one center expansion for one electron system
 */

#include <petscmat.h>
#include "fem_inf.h"
#include "angmoment.h"

struct _p_OCE1 {
  MPI_Comm comm;
  PetscReal mu;
  FEMInf fem;
  Y1s y1s;
  Mat s_r;
  Mat s_y;
};
typedef struct _p_OCE1* OCE1;

PetscErrorCode OCE1Create(OCE1 *oce1, MPI_Comm comm);
PetscErrorCode OCE1Destroy(OCE1 *oce1);
PetscErrorCode OCE1Set(OCE1 self, FEMInf fem, Y1s y1s);
PetscErrorCode OCE1CreateFromOptions(OCE1 *oce1, MPI_Comm comm);
PetscErrorCode OCE1View(OCE1 self);

PetscErrorCode OCE1GetSizes(OCE1 self, int *n_r, int *n_y);

PetscErrorCode OCE1SetSMat(OCE1 self, Mat *M);
PetscErrorCode OCE1SetSMatNullable(OCE1 self, Mat *M);
PetscErrorCode OCE1SetTMat(OCE1 self, Mat *M);
PetscErrorCode OCE1PlusPOTMat(OCE1 self, RotSym sym, POT pot, Mat M);
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat M);

#ifdef __cplusplus
}
#endif 
#endif
