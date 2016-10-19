#ifndef SCALE_H
#define SCALE_H

#include <petscmat.h>

#ifdef __cplusplus
extern "C" {
#endif 

struct _p_ScalerNone {
  MPI_Comm comm;
};
typedef struct _p_ScalerNone* ScalerNone;
PetscErrorCode ScalerNoneCreate(ScalerNone *scaler, MPI_Comm comm);
PetscErrorCode ScalerNoneDestroy(ScalerNone *scaler);
PetscErrorCode ScalerNoneView(ScalerNone self);
PetscErrorCode ScalerNoneSetRr(ScalerNone self, PetscReal xs[], int n, PetscScalar ys[]);
PetscErrorCode ScalerNoneSetQr(ScalerNone self, PetscReal xs[], int n, PetscScalar ys[]);

struct _p_ScalerUniform {
  PetscReal theta;
  MPI_Comm comm;
};
typedef struct _p_ScalerUniform* ScalerUniform;
PetscErrorCode ScalerUniformCreate(ScalerUniform *scaler, MPI_Comm comm, PetscReal theta);
PetscErrorCode ScalerUniformDestroy(ScalerUniform *scaler);
PetscErrorCode ScalerUniformView(ScalerUniform self);
PetscErrorCode ScalerUniformSetRr(ScalerUniform self, PetscReal xs[], int n, PetscScalar ys[]);
PetscErrorCode ScalerUniformSetQr(ScalerUniform self, PetscReal xs[], int n, PetscScalar ys[]);

struct _p_ScalerSharpECS {
  PetscReal r0;
  PetscReal theta;
  MPI_Comm comm;
};
typedef struct _p_ScalerSharpECS* ScalerSharpECS;
PetscErrorCode ScalerSharpECSCreate(ScalerSharpECS *scaler, MPI_Comm comm, PetscReal r0, PetscReal theta);
PetscErrorCode ScalerSharpECSDestroy(ScalerSharpECS *scaler);
PetscErrorCode ScalerSharpECSView(ScalerSharpECS self);
PetscErrorCode ScalerSharpECSSetRr(ScalerSharpECS self, PetscReal xs[], int n, PetscScalar ys[]);
PetscErrorCode ScalerSharpECSSetQr(ScalerSharpECS self, PetscReal xs[], int n, PetscScalar ys[]);

typedef struct {
  PetscErrorCode (*Destroy)();
  PetscErrorCode (*View)();
  PetscErrorCode (*SetRr)();
  PetscErrorCode (*SetQr)();
} ScalerVtbl;
struct _p_Scaler {
  ScalerVtbl *vtbl;
  void *obj;
  MPI_Comm comm;
};
typedef struct _p_Scaler* Scaler;
PetscErrorCode ScalerCreateNone(Scaler *scaler, MPI_Comm comm);
PetscErrorCode ScalerCreateUniform(Scaler *scaler, MPI_Comm comm, PetscReal theta);
PetscErrorCode ScalerCreateSharpECS(Scaler *scaler, MPI_Comm comm, PetscReal r0, PetscReal theta);
PetscErrorCode ScalerCreateFromOptions(Scaler *scaler, MPI_Comm);
PetscErrorCode ScalerDestroy(Scaler *scaler);
PetscErrorCode ScalerView(Scaler self);
PetscErrorCode ScalerSetRr(Scaler self, PetscReal xs[], int n, PetscScalar ys[]);
PetscErrorCode ScalerSetQr(Scaler self, PetscReal xs[], int n, PetscScalar ys[]);


#ifdef __cplusplus
}
#endif

#endif
