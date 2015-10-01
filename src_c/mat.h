#ifndef MAT_H
#define MAT_H

#include <stdlib.h>
#include <string.h>
#include <slepceps.h>

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat);

#endif
