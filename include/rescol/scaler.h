#ifndef SCALER_H
#define SCALER_H
#ifdef __cplusplus
extern "C" {
#endif 
/*
  arbitary complex scaling 
*/
#include <petscmat.h>
struct _p_Scaler {
  MPI_Comm comm;
  char name[20];
  int num;
  PetscScalar *vs;
  PetscErrorCode (*Calc)(struct _p_Scaler*, PetscReal, PetscScalar*, PetscScalar*);  
  PetscErrorCode (*View)(struct _p_Scaler*, PetscViewer);
};
typedef struct _p_Scaler* Scaler;

PetscErrorCode ScalerCreate(MPI_Comm comm, Scaler *p_self);
PetscErrorCode ScalerDestroy(Scaler *p_self);

PetscErrorCode ScalerView(Scaler self, PetscViewer v);
PetscErrorCode ScalerCalc(Scaler self, PetscReal* x, int n,
			  PetscScalar *qrs, PetscScalar *Rrs);

PetscErrorCode ScalerSetNone(Scaler self);
PetscErrorCode ScalerSetUniformCS(Scaler self, PetscReal theta);
PetscErrorCode ScalerSetSharpECS(Scaler self, PetscReal r0, PetscReal theta);
PetscErrorCode ScalerSetFromOptions(Scaler self);

  // PetscBool ScalerIsType(Scaler pot, char *name);


#ifdef __cplusplus
}
#endif
#endif
