#ifndef OCE_TWO_H
#define OCE_TWO_H
#ifdef __cplusplus
extern "C" {
#endif 

/*
  matrix caluclator using one center expansions for two electron systems
*/

#include <petscmat.h>
#include "fem_inf.h"
#include "angmoment.h"

struct _p_OCE2 {

  // const
  MPI_Comm comm;
  
  // calculator
  FEMInf fem;
  Y2s y2s;

  // buffer
  Mat s_r1;
  Mat s_y2;  
};
typedef struct _p_OCE2* OCE2;

PetscErrorCode OCE2Create(OCE2 *oce2, MPI_Comm comm);
PetscErrorCode OCE2Destroy(OCE2 *oce2);
PetscErrorCode OCE2Set(OCE2 self, FEMInf fem, Y2s y2s);
PetscErrorCode OCE2CreateFromOptions(OCE2 *oce2, MPI_Comm comm);
PetscErrorCode OCE2View(OCE2 self);

PetscErrorCode OCE2SetSMat(OCE2 self, Mat *M);
PetscErrorCode OCE2SetTMat(OCE2 self, Mat *M);
PetscErrorCode OCE2PlusVneMat(OCE2 self, PetscReal a, PetscReal z, Mat *M);
PetscErrorCode OCE2PlusVeeMat(OCE2 self, Mat *M);

#ifdef __cplusplus
}
#endif 
#endif
