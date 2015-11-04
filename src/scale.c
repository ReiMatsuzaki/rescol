#include "scale.h"

PetscErrorCode ScalerNoneCreate(ScalerNone *scaler, MPI_Comm comm) {

  ScalerNone _scaler;
  PetscNew(&_scaler);

  _scaler->comm = comm;

  *scaler = _scaler;
  return 0;

}
PetscErrorCode ScalerNoneDestroy(ScalerNone *scaler) {
  PetscFree(*scaler);
  return 0;
}
PetscErrorCode ScalerNoneView(ScalerNone self) {
  PetscPrintf(self->comm, ">>>> ScalerNone >>>>\n");
  PetscPrintf(self->comm, "<<<< ScalerNone <<<<\n");
  return 0;
}
PetscErrorCode ScalerNoneSetRr(ScalerNone self, PetscReal xs[], int n, PetscScalar ys[]) {

  for(int i = 0; i < n; i++)
    ys[i] = xs[i];
  return 0;
}
PetscErrorCode ScalerNoneSetQr(ScalerNone self, PetscReal xs[], int n, PetscScalar ys[]) {

  for(int i = 0; i < n; i++)
    ys[i] = 1.0;
  return 0;

}

PetscErrorCode ScalerUniformCreate(ScalerUniform *scaler, MPI_Comm comm, PetscReal theta) {

  ScalerUniform _scaler;
  PetscNew(&_scaler);

  _scaler->theta = theta;
  _scaler->comm = comm;

  *scaler = _scaler;
  return 0;

}
PetscErrorCode ScalerUniformDestroy(ScalerUniform *scaler) {
  PetscFree(*scaler);
  return 0;
}
PetscErrorCode ScalerUniformView(ScalerUniform self) {
  PetscPrintf(self->comm, ">>>> ScalerUniform >>>>\n");
  PetscPrintf(self->comm, "theta: %f\n", self->theta);
  PetscPrintf(self->comm, "<<<< ScalerUniform <<<<\n");
  return 0;
}
PetscErrorCode ScalerUniformSetRr(ScalerUniform self, PetscReal xs[], int n, PetscScalar ys[]) {

  for(int i = 0; i < n; i++)
    ys[i] = xs[i] * cos(-self->theta) + xs[i] * sin(-self->theta) * PETSC_i;
  return 0;
}
PetscErrorCode ScalerUniformSetQr(ScalerUniform self, PetscReal xs[], int n, PetscScalar ys[]) {

  for(int i = 0; i < n; i++)
    ys[i] = cos(-self->theta) + sin(-self->theta) * PETSC_i;
  return 0;

}

PetscErrorCode ScalerSharpECSCreate(ScalerSharpECS *scaler, MPI_Comm comm, 
				    PetscReal r0, PetscReal theta) {

  ScalerSharpECS _scaler;
  PetscNew(&_scaler);

  _scaler->theta = theta;
  _scaler->r0 = r0;
  _scaler->comm = comm;

  *scaler = _scaler;
  return 0;
}
PetscErrorCode ScalerSharpECSDestroy(ScalerSharpECS *scaler) {
  PetscFree(*scaler);
  return 0;
}
PetscErrorCode ScalerSharpECSView(ScalerSharpECS self) {
  PetscPrintf(self->comm, ">>>> ScalerSharpECS >>>>\n");
  PetscPrintf(self->comm, "theta: %f\n", self->theta);
  PetscPrintf(self->comm, "<<<< ScalerSharpECS <<<<\n");
  return 0;
}
PetscErrorCode ScalerSharpECSSetRr(ScalerSharpECS self, PetscReal xs[], int n, PetscScalar ys[]) {

  PetscReal r0 = self->r0;
  PetscReal t = self->theta;
  for(int i = 0; i < n; i++) {
    
    PetscReal x = xs[i];
    if(x < r0)
      ys[i] = x;
    else
      ys[i] = r0 + (x-r0) * cos(-t) + (x-r0) * sin(-t) * PETSC_i;
  }
  return 0;
}
PetscErrorCode ScalerSharpECSSetQr(ScalerSharpECS self, PetscReal xs[], int n, PetscScalar ys[]) {

  PetscReal r0 = self->r0;
  PetscReal t = self->theta;
  for(int i = 0; i < n; i++) {
    PetscReal x = xs[i];
    if(x < r0)
      ys[i] = 1.0;
    else
      ys[i] = cos(-t) + sin(-t) * PETSC_i;
  }
  return 0;

}

PetscErrorCode ScalerCreate(Scaler *scaler, MPI_Comm comm) {

  Scaler _scaler;
  PetscNew(&_scaler);
  
  _scaler->vtbl = NULL;
  _scaler->obj = NULL;
  _scaler->comm = comm;

  *scaler = _scaler;
  return 0;
}
ScalerVtbl ScalerNone_vtbl;
PetscErrorCode ScalerCreateNone(Scaler *scaler, MPI_Comm comm) {

  ScalerCreate(scaler, comm);
  ScalerNone none; ScalerNoneCreate(&none, comm);
    
  static int init = 0;
  if(init == 0) {
    ScalerNone_vtbl.Destroy =ScalerNoneDestroy;
    ScalerNone_vtbl.View =   ScalerNoneView;
    ScalerNone_vtbl.SetRr =  ScalerNoneSetRr;
    ScalerNone_vtbl.SetQr =  ScalerNoneSetQr;
  }

  (*scaler)->vtbl = &ScalerNone_vtbl;
  (*scaler)->obj = none;

  return 0;
}
ScalerVtbl ScalerUniform_vtbl;
PetscErrorCode ScalerCreateUniform(Scaler *scaler, MPI_Comm comm, PetscReal theta) {
  
  ScalerCreate(scaler, comm);
  ScalerUniform uni; ScalerUniformCreate(&uni, comm, theta);
  
  static int init = 0;
  if(init == 0) {
    ScalerUniform_vtbl.Destroy =ScalerUniformDestroy;
    ScalerUniform_vtbl.View = ScalerUniformView;
    ScalerUniform_vtbl.SetRr = ScalerUniformSetRr;
    ScalerUniform_vtbl.SetQr = ScalerUniformSetQr;
    
  }

  (*scaler)->vtbl = &ScalerUniform_vtbl;
  (*scaler)->obj = uni;

  return 0;
}
ScalerVtbl ScalerSharpECS_vtbl;
PetscErrorCode ScalerCreateSharpECS(Scaler *scaler, MPI_Comm comm, PetscReal r0, PetscReal theta) {
  
  ScalerCreate(scaler, comm);
  ScalerSharpECS ecs; 
  ScalerSharpECSCreate(&ecs, comm, r0, theta);
  
  static int init = 0;
  if(init == 0) {
    ScalerSharpECS_vtbl.Destroy =ScalerSharpECSDestroy;
    ScalerSharpECS_vtbl.View = ScalerSharpECSView;
    ScalerSharpECS_vtbl.SetRr = ScalerSharpECSSetRr;
    ScalerSharpECS_vtbl.SetQr = ScalerSharpECSSetQr;
    
  }

  (*scaler)->vtbl = &ScalerSharpECS_vtbl;
  (*scaler)->obj = ecs;

  return 0;
}
PetscErrorCode ScalerCreateFromOptions(Scaler *scaler, MPI_Comm comm) {

  PetscErrorCode ierr;
  char type[10] = "none";
  PetscReal r0 = 0.0;
  PetscReal theta = 0.0;

  ierr = PetscOptionsGetString(NULL, "-scaler_type", type, 10, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, "-scaler_r0", &r0, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, "-scaler_theta", &theta, NULL); CHKERRQ(ierr);

  if(strcmp(type, "none") == 0) {
    ScalerCreateNone(scaler, comm);
  } else if(strcmp(type, "uni") == 0) {
    ScalerCreateUniform(scaler, comm, theta*M_PI/180.0);
  } else if(strcmp(type, "secs") == 0) {
    ScalerCreateSharpECS(scaler, comm, r0, theta);
  } else
    SETERRQ(comm, 1, "scaler_type<-{none, uni, secs}");
  
  return 0;
}
PetscErrorCode ScalerDestroy(Scaler *scaler) {
  Scaler self = *scaler;
  if(self->vtbl->Destroy == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->vtbl->Destroy(&self->obj);
  return 0;
}
PetscErrorCode ScalerView(Scaler self) {
  if(self->vtbl->View == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->vtbl->View(self->obj);
  return 0;
}
PetscErrorCode ScalerSetRr(Scaler self, PetscReal xs[], int n, PetscScalar ys[]) {
  if(self->vtbl->SetRr == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->vtbl->SetRr(self->obj, xs, n, ys);
  return 0;
}
PetscErrorCode ScalerSetQr(Scaler self, PetscReal xs[], int n, PetscScalar ys[]) {
  if(self->vtbl->SetQr == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->vtbl->SetQr(self->obj, xs, n, ys);
  return 0;
}
