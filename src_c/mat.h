#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <string.h>
#include <slepceps.h>

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);
PetscErrorCode EPSWriteToFile(EPS eps, char* path_detail, char* path_eigvals, char* path_eigvecs);

#endif
