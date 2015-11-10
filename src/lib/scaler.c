#include <rescol/scaler.h>

PetscErrorCode ScalerCreate(MPI_Comm comm, Scaler *p_self) {
  Scaler self;
  PetscNew(&self);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode ScalerDestroy(Scaler *p_self) {
  PetscErrorCode ierr;
  ierr = PetscFree((*p_self)->vs); CHKERRQ(ierr);
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode ScalerView(Scaler self, PetscViewer v) {
  PetscErrorCode ierr;

  PetscViewerType type;
  PetscViewerGetType(v, &type);
  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "only ascii is supported");

  PetscViewerASCIIPrintf(v, ">>>> Scaler >>>>\n");
  PetscViewerASCIIPrintf(v, "name: %s\n", self->name);  
  ierr = self->View(self, v);
  PetscViewerASCIIPrintf(v, "<<<< Scaler <<<<\n");
  
  return 0;
}
PetscErrorCode ScalerCalc(Scaler self, PetscReal* xs, int n, 
			  PetscScalar* qrs, PetscScalar* Rrs) {
  for(int i = 0; i < n; i++) {
    PetscReal x = xs[i];
    PetscScalar q, R;
    self->Calc(self, x, &q, &R);
    if(qrs != NULL)
      qrs[i] = q;
    if(Rrs != NULL)
      Rrs[i] = R;
  }
  return 0;
}

PetscErrorCode SharpECSView(Scaler self, PetscViewer v) {

  PetscReal r0 = self->vs[0];
  PetscReal t  = self->vs[1];

  PetscViewerASCIIPrintf(v, "q(x) = | 1      (x < r0) \n");
  PetscViewerASCIIPrintf(v, "       | e^(it) (x > r0) \n");
  PetscViewerASCIIPrintf(v, "r0: %f\n", r0);
  PetscViewerASCIIPrintf(v, "t: %f\n", t);

  return 0;
}
PetscErrorCode SharpECSCacl(Scaler self, PetscReal x, 
			    PetscScalar *q, PetscScalar *R) {

  PetscReal r0 = self->vs[0];
  PetscReal t  = self->vs[1];

  if(x < r0) {    
    *q = 1.0;
    *R = x;
  } else {
    PetscScalar expi = cos(t) + sin(t)*PETSC_i;
    *q = expi;
    *R = r0 + (x-r0) * expi;
  }

  return 0;
}

PetscErrorCode ScalerSetNone(Scaler self) {
  PetscErrorCode ierr;
  ierr = ScalerSetUniformCS(self, 0.0); CHKERRQ(ierr);
  strcpy(self->name, "None");
  return 0;
}
PetscErrorCode ScalerSetUniformCS(Scaler self, PetscReal t) {
  PetscErrorCode ierr;
  ierr = ScalerSetSharpECS(self, 0.0, t); CHKERRQ(ierr);
  strcpy(self->name, "UniformCS");
  return 0;
}
PetscErrorCode ScalerSetSharpECS(Scaler self, PetscReal r0, PetscReal t) {

  strcpy(self->name, "SharpECS");
  self->num = 2;
  PetscMalloc1(2, &self->vs);
  self->vs[0] = r0;
  self->vs[1] = t;
  self->Calc = SharpECSCacl;
  self->View = SharpECSView;
  return 0;
}
PetscErrorCode ScalerSetFromOptions(Scaler self) {

  PetscErrorCode ierr;
  char type[10] = "none";

  ierr = PetscOptionsGetString(NULL, "-scaler_type", type, 10, NULL); CHKERRQ(ierr);

  if(strcmp(type, "none") == 0) {
    ScalerSetNone(self);
  } else if(strcmp(type, "uniformCS") == 0) {
    PetscReal t = 0.0;
    ierr = PetscOptionsGetReal(NULL, "-scaler_theta", &t, NULL); CHKERRQ(ierr);
    ScalerSetUniformCS(self, t*M_PI/180.0);
  } else if(strcmp(type, "sharpECS") == 0) {
    PetscReal r0 = 0.0;
    PetscReal t = 0.0;
    ierr = PetscOptionsGetReal(NULL, "-scaler_r0", &r0, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, "-scaler_theta", &t, NULL); CHKERRQ(ierr);
    ScalerSetSharpECS(self, r0, t*M_PI/180.0);
  } else
    SETERRQ(self->comm, 1, "scaler_type<-{none, uniformCS, sharpECS}");
  return 0;  
}


