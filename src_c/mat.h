#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <string.h>
#include <slepceps.h>

PetscErrorCode MatCreateFromCOOFormatFileHandler(FILE* path, Mat* mat);
PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);
PetscErrorCode EPSWriteToFile(EPS eps, char* path_detail, char* path_eigvals, char* path_eigvecs);

PetscErrorCode MatInitSynthesize(Mat A, Mat B, MPI_Comm comm, Mat *C);
PetscErrorCode MatSynthesize(Mat A, Mat B, PetscScalar c, Mat *C, InsertMode mode);
PetscErrorCode MatSetSynthesize(Mat A, Mat B, PetscScalar c, MPI_Comm comm, Mat *C);

PetscErrorCode MatInitSynthesize3(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D);
PetscErrorCode MatSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, Mat *D, InsertMode mode);
PetscErrorCode MatSetSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, MPI_Comm comm, Mat *D);

#endif
