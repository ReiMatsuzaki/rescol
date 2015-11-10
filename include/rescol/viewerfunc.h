#ifndef FUNCVIEWER_H
#define FUNCVIEWER_H
#ifdef __cplusplus
// extern "C" {
#endif 
#include <petscviewer.h>

struct _p_ViewerFunc {
  PetscViewer base;
  MPI_Comm comm;
  int num;
  PetscReal *xs;
  PetscBool active_base;
  PetscBool active_range;
};
typedef struct _p_ViewerFunc* ViewerFunc;
PetscErrorCode ViewerFuncCreate(MPI_Comm comm, ViewerFunc *p_self);
PetscErrorCode ViewerFuncDestroy(ViewerFunc *p_self);

PetscErrorCode ViewerFuncView(ViewerFunc self, PetscViewer viewer);

PetscErrorCode ViewerFuncSetBase(ViewerFunc self, PetscViewer base);
PetscErrorCode ViewerFuncSetRange(ViewerFunc self, int num, PetscReal xmax);
PetscErrorCode ViewerFuncSetFromOptions(ViewerFunc self);

PetscErrorCode ViewerFuncGetXs(ViewerFunc self, int *num, PetscReal **xs);
PetscBool ViewerFuncIsActive(ViewerFunc self);
PetscErrorCode ViewerFuncCheckAcrive(ViewerFunc self);

#ifdef __cplusplus
  // }
#endif
#endif
