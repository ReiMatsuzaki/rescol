#ifndef FUNCVIEWER_H
#define FUNCVIEWER_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <petscviewer.h>

  /*
    Options data base
    .  -{prefix}viewerfunc_path  : path to write function values
    .  -{prefix}viewerfunc_num   : number of grid
    .  -{prefix}viewerfunc_xmax  : maximum of grid
   */
  
struct _p_ViewerFunc {
  PetscViewer base;
  MPI_Comm comm;
  int num;
  PetscReal *xs;

  PetscBool active_base;
  PetscBool active_range;

  char opt_prefix[100];
  char opt_base[100];
};
typedef struct _p_ViewerFunc* ViewerFunc;
PetscErrorCode ViewerFuncCreate(MPI_Comm comm, ViewerFunc *p_self);
PetscErrorCode ViewerFuncDestroy(ViewerFunc *p_self);

PetscViewer ViewerFuncGetBase();
PetscErrorCode ViewerFuncView(ViewerFunc self, PetscViewer viewer);

PetscErrorCode ViewerFuncSetBase(ViewerFunc self, PetscViewer base);
PetscErrorCode ViewerFuncSetRange(ViewerFunc self, int num, PetscReal xmin, PetscReal xmax);
PetscErrorCode ViewerFuncSetOptionsPrefix(ViewerFunc self, const char prefix[]);
PetscErrorCode ViewerFuncSetFromOptions(ViewerFunc self);

PetscErrorCode ViewerFuncGetRange(ViewerFunc self, int *num, PetscReal **xs);
PetscBool ViewerFuncIsActive(ViewerFunc self);

#ifdef __cplusplus
 }
#endif
#endif
