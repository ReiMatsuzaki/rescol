#ifndef SYNTHESIZE_H
#define SYNTHESIZE_H
#ifdef __cplusplus
extern "C" {
#endif 

#include <petscmat.h>
  //#include <src/mat/impls/aij/seq/aij.h> /*I "petscmat.h" I*/

PetscErrorCode VecVecSynthesizeSymbolic(Vec A, Vec B, Vec *C);
PetscErrorCode VecVecSynthesizeNumeric(Vec A, Vec B, PetscScalar c, Vec C);
PetscErrorCode VecVecSynthesize(Vec A, Vec B, PetscScalar a, MatReuse scall, Vec *C);

PetscErrorCode MatMatSynthesizeSymbolic(Mat A, Mat B, Mat *C);
PetscErrorCode MatMatSynthesizeNumeric(Mat A, Mat B, PetscScalar c, Mat C);
PetscErrorCode MatMatSynthesize(Mat A, Mat B, PetscScalar a, MatReuse scall, Mat *C);

PetscErrorCode MatMatMatSynthesizeSymbolic(Mat A, Mat B, Mat C, Mat *D);
PetscErrorCode MatMatMatSynthesizeNumeric(Mat A, Mat B, Mat C, PetscScalar d, Mat D);
PetscErrorCode MatMatMatSynthesize(Mat A, Mat B, Mat C, PetscScalar a, MatReuse scall, Mat *D);

#ifdef __cplusplus
}
#endif
#endif
