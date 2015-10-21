#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <slepceps.h>

PetscErrorCode MatCreateFromCOOFormatFileHandler(FILE* path, Mat* mat);
PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);

PetscErrorCode MatSetDirFile(const char* dn, const char* fn, Mat *M);
PetscErrorCode VecCreateFromFile(const char* path, MPI_Comm comm, Vec *v );
PetscErrorCode PrintTimeStamp(MPI_Comm comm, const char* label, time_t *t);
PetscErrorCode EPSWriteToFile(EPS eps, char* path_detail, char* path_eigvals, char* path_eigvecs);

PetscErrorCode VecInitSynthesize(Vec A, Vec B, MPI_Comm comm, Vec *C);
PetscErrorCode VecSynthesize(Vec A, Vec B, PetscScalar c, Vec *C, InsertMode mode);
PetscErrorCode VecSetSynthesize(Vec A, Vec B, PetscScalar c, MPI_Comm comm, Vec *C);

PetscErrorCode MatInitSynthesize(Mat A, Mat B, MPI_Comm comm, Mat *C);
PetscErrorCode MatSynthesize(Mat A, Mat B, PetscScalar c, Mat *C, InsertMode mode);
PetscErrorCode MatSetSynthesize(Mat A, Mat B, PetscScalar c, MPI_Comm comm, Mat *C);
PetscErrorCode MatSetSynthesizeSlow(Mat A, Mat B, PetscScalar c, MPI_Comm comm, Mat *C);
PetscErrorCode MatSetSynthesizeFast(Mat A, Mat B, MPI_Comm comm, Mat *C);

PetscErrorCode MatInitSynthesize3(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D);
PetscErrorCode MatSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, Mat *D, InsertMode mode);
PetscErrorCode MatSetSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, MPI_Comm comm, Mat *D);
PetscErrorCode MatSetSynthesize3Fast(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D);

#endif
