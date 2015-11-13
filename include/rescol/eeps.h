#ifndef EEPS_H
#define EEPS_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <slepceps.h>

/*
  enhanced EPS object. Support additional these properties:
  1. Read Guess vector option
  2. Write more detail for eigen values and its error
*/

struct _p_EEPS {
  EPS eps;
  Mat S;
  PetscViewer viewer_values;
};
typedef struct _p_EEPS* EEPS;
PetscErrorCode EEPSCreate(MPI_Comm comm, EEPS *p_self);
PetscErrorCode EEPSDestroy(EEPS *p_self);

PetscErrorCode EEPSSetOperators(EEPS self, Mat H, Mat S );
PetscErrorCode EEPSSetTarget(EEPS self, PetscScalar target);
PetscErrorCode EEPSSetInitSpace(EEPS self, int num_guess, Vec *guess);
PetscErrorCode EEPSSetFromOptions(EEPS self);

PetscErrorCode EEPSSolve(EEPS self);
PetscErrorCode EEPSGetEigenvector(EEPS self, int i, Vec k);

#ifdef __cplusplus
}
#endif
#endif
