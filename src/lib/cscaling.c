#include <petscviewer.h>
#include <rescol/cscaling.h>


PetscErrorCode NoneApply(void *ctx, PetscInt n,
			 const PetscScalar *xs, PetscScalar *ys) {
  ys[0] = 1.0;
  ys[1] = xs[0];
  return 0;
}
PetscErrorCode NoneView(void *ctx, PetscViewer v) {
  PetscViewerASCIIPrintf(v, "Complex Scaling (None)\n");
  return 0;
}

typedef struct {
  PetscReal theta;
} UniformCS;
PetscErrorCode UniformCSApply(void *ctx, PetscInt n,
			      const PetscScalar *xs, PetscScalar *ys) {
  UniformCS* cs = (UniformCS*)ctx;
  PetscReal t = cs->theta;  
  ys[0] = cos(t) + PETSC_i * sin(t);
  ys[1] = xs[0] * ys[0];
  return 0;
}
PetscErrorCode UniformCSView(void *ctx, PetscViewer v) {
  UniformCS* cs = (UniformCS*)ctx;
  PetscReal t = cs->theta;    
  PetscViewerASCIIPrintf(v, "Complex Scaling (Uniform)\n");
  PetscViewerASCIIPrintf(v, "theta : %f\n", t);
  return 0;
}

typedef struct {
  PetscReal r0;
  PetscReal theta;
} SharpECS;
PetscErrorCode SharpECSApply(void *ctx, PetscInt n,
			      const PetscScalar *xs, PetscScalar *ys) {
  SharpECS* cs = (SharpECS*)ctx;
  PetscReal r0 = cs->r0;  
  PetscReal t = cs->theta;  
  PetscReal x  = xs[0];
  if(x < r0) {
    ys[0] = 1.0;
    ys[1] = x;
  } else {
    ys[0] = cos(t) + PETSC_i * sin(t);
    ys[1] = r0 + (xs[0]-r0)*ys[0];
  }
  return 0;
}
PetscErrorCode SharpECSView(void *ctx, PetscViewer v) {
  SharpECS* cs = (SharpECS*)ctx;
  PetscReal r0 = cs->r0;  
  PetscReal t = cs->theta;  
  PetscViewerASCIIPrintf(v, "Complex Scaling (Sharp Exterior)\n");
  PetscViewerASCIIPrintf(v, "r0    = %f\n", r0);
  PetscViewerASCIIPrintf(v, "theta = %f\n", t);
  return 0;
}

PetscErrorCode CScalingCreate(MPI_Comm comm, CScaling *p_self) {
  PFCreate(comm, 1, 2, p_self);
  return 0;
}
PetscErrorCode CScalingSetNone(CScaling self) {
  PFSet(self, NoneApply, NULL, NoneView, NULL, NULL);
  return 0;
}
PetscErrorCode CScalingSetUniformCS(CScaling self, PetscReal t) {
  UniformCS *ctx; PetscNew(&ctx);
  ctx->theta = t;
  PFSet(self, UniformCSApply, NULL, UniformCSView, NULL, ctx);
  return 0;  
}
PetscErrorCode CScalingSetSharpECS(CScaling self, PetscReal r0, PetscReal t) {
  SharpECS *ctx; PetscNew(&ctx);
  ctx->r0 = r0;
  ctx->theta = t;
  PFSet(self, SharpECSApply, NULL, SharpECSView, NULL, ctx);
  return 0;    
}
PetscErrorCode CScalingSetFromOptions(CScaling self) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  char type[10] = "none";
  PetscReal r0 = 0.0;
  PetscReal theta = 0.0;

  PetscOptionsGetString(NULL, "-cscaling_type", type, 10, NULL);
  PetscOptionsGetReal(NULL, "-cscaling_r0", &r0, NULL);
  PetscOptionsGetReal(NULL, "-cscaling_theta", &theta, NULL); 

  if(strcmp(type, "none") == 0) {
    CScalingSetNone(self);
  } else if(strcmp(type, "uniform_cs") == 0) {
    CScalingSetUniformCS(self, theta*M_PI/180.0);
  } else if(strcmp(type, "sharp_ecs") == 0) {
    CScalingSetSharpECS(self, r0, theta*M_PI/180.0);
  } else
    SETERRQ(comm, 1, "cscaling_type<-{none, uni, secs}");
  return 0;
}

PetscErrorCode CScalingCalc(CScaling self, PetscReal *xs, int n,
			    PetscScalar *qrs, PetscScalar *Rrs) {
  for(int i = 0; i < n; i++) {
    PetscScalar x[1] = {xs[i]};
    PetscScalar ys[2];
    PFApply(self, 1, x, ys);
    if(qrs!=NULL)
      qrs[i] = ys[0];
    if(Rrs!=NULL)
      Rrs[i] = ys[1];
  }
  return 0;
}


