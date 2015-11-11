#ifndef MAT_H
#define MAT_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <slepceps.h>

PetscReal ScalarAbs(PetscScalar x);

PetscErrorCode MatCreateFromCOOFormatFileHandler(FILE* path, Mat* mat);
PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);

PetscErrorCode MatSetDirFile(const char* dn, const char* fn, Mat *M);
PetscErrorCode VecCreateFromFile(const char* path, MPI_Comm comm, Vec *v );
PetscErrorCode PrintTimeStamp(MPI_Comm comm, const char* label, time_t *t);
PetscErrorCode EPSWriteToFile(EPS eps, char* path_detail, char* path_eigvals, char* path_eigvecs);
PetscErrorCode EPSCreateForBoundState(EPS *eps, MPI_Comm comm, Mat H, Mat S, PetscScalar target);
PetscErrorCode EPSSetGuessFromFiles(EPS eps, MPI_Comm comm, char **fn_list, int n);

PetscErrorCode VecSynthesizeSymbolic(Vec A, Vec B, Vec *C);
PetscErrorCode VecSynthesizeNumeric(Vec A, Vec B, PetscScalar c, Vec C);
PetscErrorCode VecSynthesize(Vec A, Vec B, PetscScalar c, MatReuse scall, Vec *C);
PetscErrorCode VecNormalizeForS(Mat S, Vec x);

PetscErrorCode MatSynthesizeSymbolic(Mat A, Mat B, Mat *C);
PetscErrorCode MatSynthesizeNumeric(Mat A, Mat B, PetscScalar c, Mat C);
PetscErrorCode MatSynthesize(Mat A, Mat B, PetscScalar c, MatReuse scall, Mat *C);
  
  //PetscErrorCode MatInitSynthesize(Mat A, Mat B, MPI_Comm comm, Mat *C);
  //PetscErrorCode MatSynthesizeOld(Mat A, Mat B, PetscScalar c, Mat *C, InsertMode mode);
  //PetscErrorCode MatSetSynthesize(Mat A, Mat B, PetscScalar c, MPI_Comm comm, Mat *C);
  //PetscErrorCode MatSetSynthesizeSlow(Mat A, Mat B, PetscScalar c, MPI_Comm comm, Mat *C);
  //PetscErrorCode MatSetSynthesizeFast(Mat A, Mat B, MPI_Comm comm, Mat *C);

PetscErrorCode MatSynthesize3Symbolic(Mat A, Mat B, Mat C, Mat *D);
PetscErrorCode MatSynthesize3Numeric(Mat A, Mat B, Mat C, PetscScalar d, Mat D);
PetscErrorCode MatSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, MatReuse scall, Mat *D);
  //PetscErrorCode MatInitSynthesize3(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D);
  //PetscErrorCode MatSynthesize3Old(Mat A, Mat B, Mat C, PetscScalar d, Mat *D, InsertMode mode);
  PetscErrorCode MatSetSynthesize3Old(Mat A, Mat B, Mat C, PetscScalar d, MPI_Comm comm, Mat *D);
  //PetscErrorCode MatSetSynthesize3Fast(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D);

PetscErrorCode PartialCoulomb(int q, double r1, double r2, double *y);
PetscErrorCode LegGauss(int n, int i, PetscScalar *x, PetscScalar *w);
PetscErrorCode LobGauss(int n, int i, PetscScalar *x, PetscScalar *w);

#ifdef __cplusplus
}
#endif
#endif
