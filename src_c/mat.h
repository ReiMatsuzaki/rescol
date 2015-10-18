#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <string.h>
#include <slepceps.h>

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);
PetscErrorCode EPSWriteToFile(EPS eps, const char* path_detail, const char* path_eigvals, const char* path_eigvecs);

#endif
