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
#include "y2s.h"
#include "synthesize.h"

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

PetscErrorCode OCE2Create(MPI_Comm comm, OCE2 *p_self);
PetscErrorCode OCE2Destroy(OCE2 *p_self);

PetscErrorCode OCE2View(OCE2 self, PetscViewer v);

PetscErrorCode OCE2Set(OCE2 self, FEMInf fem, Y2s y2s);
PetscErrorCode OCE2SetFromOptions(OCE2 self);

PetscErrorCode OCE2GetSizes(OCE2 self, int *n_r1, int *n_y2);

PetscErrorCode OCE2CreateMat(OCE2 self, Mat *M);
PetscErrorCode OCE2SMat(OCE2 self, MatReuse scall, Mat *M, PetscBool *is_id);
PetscErrorCode OCE2TMat(OCE2 self, MatReuse scall, Mat *M);
PetscErrorCode OCE2VneMat(OCE2 self, PetscReal a, PetscReal z, 
			  MatReuse scall, Mat *M);
PetscErrorCode OCE2VeeMat(OCE2 self, MatReuse scall, Mat *M);
PetscErrorCode OCE2PlusVneMat(OCE2 self, PetscReal a, PetscReal z, Mat M);
PetscErrorCode OCE2PlusVeeMat(OCE2 self, Mat M);
#ifdef __cplusplus
}
#endif 
#endif
