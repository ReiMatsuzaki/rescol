#include <petscviewer.h>
#include "../include/math.h"
#include "../include/cscaling.h"


PetscErrorCode NoneApply(void *ctx, PetscInt n,
			 const PetscScalar *xs, PetscScalar *ys) {
  ys[0] = 1.0;
  ys[1] = xs[0];
  return 0;
}
PetscErrorCode NoneView(void *ctx, PetscViewer v) {
  PetscViewerASCIIPrintf(v, "type: Complex Scaling (None)\n");
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
  PetscViewerASCIIPrintf(v, "type :  Complex Scaling (Uniform)\n");
  PetscViewerASCIIPrintf(v, "theta : %f\n", t);
  return 0;
}
PetscErrorCode UniformCSDestroy(void *ctx) {
  UniformCS* cs = (UniformCS*)ctx;  
  PetscFree(cs);
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
  PetscViewerASCIIPrintf(v, "type  : Complex Scaling (Sharp Exterior)\n");
  PetscViewerASCIIPrintf(v, "r0    : %f\n", r0);
  PetscViewerASCIIPrintf(v, "theta : %f\n", t);
  return 0;
}
PetscErrorCode SharpECSDestroy(void *ctx) {
  SharpECS* cs = (SharpECS*)ctx;
  PetscFree(cs);
  return 0;
}

// ---- Basics ----
PetscErrorCode CScalingCreate(MPI_Comm comm, CScaling *p_self) {
  
  CScaling self;
  PetscMalloc1(1, &self);
  self->comm = comm;
  self->use_cscaling = PETSC_FALSE;
  PFCreate(comm, 1, 2, &self->pf);
  self->R0 = 0.0;

  *p_self = self;
  return 0;
}
PetscErrorCode CScalingCopy(CScaling self, CScaling other) {
  PetscErrorCode ierr;
  if(self->type == CScalingNone) {
    ierr = CScalingSetNone(other); CHKERRQ(ierr);
  } else if(self->type == CScalingUniformCS) {
    PetscScalar t = ((UniformCS*)self->ctx)->theta;
    ierr = CScalingSetUniformCS(other, t); CHKERRQ(ierr);
  } else if(self->type == CScalingSharpECS) {
    SharpECS* ctx = (SharpECS*)self->ctx;
    ierr = CScalingSetSharpECS(other, ctx->r0, ctx->theta); CHKERRQ(ierr);
  } else {
    SETERRQ(self->comm, 1, "Invalid type");
  }
  return 0;

}
PetscErrorCode CScalingDestroy(CScaling *p_self) {
  CScaling self = *p_self;
  PFDestroy(&self->pf);
  PetscFree(*p_self);
  return 0;
}
PetscErrorCode CScalingView(CScaling self, PetscViewer v) {
  
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  PetscViewerType type;     PetscViewerGetType(v, &type);
  PetscViewerFormat format; PetscViewerGetFormat(v, &format);

  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);    

  if(iascii) {
    PetscViewerASCIIPrintf(v, "CScaling object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "use_cscaling: %s\n",
			   (self->use_cscaling? "Yes" : "No"));
    if(self->use_cscaling) {
      PetscViewerASCIIPrintf(v, "R0:    %f\n", self->R0);
      PetscViewerASCIIPrintf(v, "theta: %f\n", self->theta);
      PetscViewerASCIIPrintf(v, "pf:\n");
      PFView(self->pf, v);
    }
    PetscViewerASCIIPopTab(v);
  }
  return 0;
}
PetscErrorCode CScalingSetNone(CScaling self) {
  PFSet(self->pf, NoneApply, NULL, NoneView, NULL, NULL);
  self->use_cscaling = PETSC_FALSE;
  self->R0 = 0.0;
  self->theta = 0.0;
  self->type = CScalingNone;
  self->ctx = NULL;
  return 0;
}
PetscErrorCode CScalingSetUniformCS(CScaling self, PetscReal t) {
  UniformCS *ctx; PetscNew(&ctx);
  ctx->theta = t;
  PFSet(self->pf, UniformCSApply, NULL, UniformCSView, UniformCSDestroy, ctx);
  self->use_cscaling = PETSC_TRUE;
  self->R0 = 0.0;
  self->theta = t;
  self->type = CScalingUniformCS;
  self->ctx = ctx;
  return 0;  
}
PetscErrorCode CScalingSetSharpECS(CScaling self, PetscReal r0, PetscReal t) {
  SharpECS *ctx; PetscNew(&ctx);
  ctx->r0 = r0;
  ctx->theta = t;
  PFSet(self->pf, SharpECSApply, NULL, SharpECSView, SharpECSDestroy, ctx);
  self->use_cscaling = PETSC_TRUE;
  self->R0 = r0;
  self->theta = t;
  self->type = CScalingSharpECS;
  self->ctx = ctx;
  return 0;    
}
PetscErrorCode CScalingSetFromOptions(CScaling self) {

  PetscErrorCode ierr;
  char type[10] = "none";
  PetscReal r0 = 0.0;
  PetscReal theta = 0.0;

  ierr = PetscOptionsGetString(NULL,NULL,"-cscaling_type", type, 10, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-cscaling_r0", &r0, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-cscaling_theta", &theta, NULL);  CHKERRQ(ierr);

  if(strcmp(type, "none") == 0) {
    CScalingSetNone(self);
  } else if(strcmp(type, "uniform_cs") == 0) {
    CScalingSetUniformCS(self, theta*M_PI/180.0);
  } else if(strcmp(type, "sharp_ecs") == 0) {
    CScalingSetSharpECS(self, r0, theta*M_PI/180.0);
  } else
    SETERRQ(self->comm, 1, "cscaling_type<-{none, uniform_cs, sharp_ecs}");
  return 0;
}

// ---- Calculation ----
PetscErrorCode CScalingCalc(CScaling self, PetscReal *xs, int n,
			    PetscScalar *qrs, PetscScalar *Rrs) {
  for(int i = 0; i < n; i++) {
    PetscScalar x[1] = {xs[i]};
    PetscScalar ys[2];
    PFApply(self->pf, 1, x, ys);
    if(qrs!=NULL)
      qrs[i] = ys[0];
    if(Rrs!=NULL)
      Rrs[i] = ys[1];
  }
  return 0;
}
PetscErrorCode CScalingCalcOne(CScaling self, PetscReal x,
			       PetscScalar *qr, PetscScalar *Rr) {
  PetscScalar in_x[1] = {x};
  PetscScalar out_y[2];
  PFApply(self->pf, 1, in_x, out_y);
  if(qr != NULL)
    *qr = out_y[0];
  if(Rr != NULL)
    *Rr = out_y[1];
  return 0;
}
PetscErrorCode CScalingQ(CScaling self, PetscBool *_use_cscaling) {
  *_use_cscaling = self->use_cscaling;
  return 0;
}
PetscErrorCode CScalingGetRadius(CScaling self, PetscReal *R0) {
  *R0 = self->R0;
  return 0;
}
PetscErrorCode CScalingGetTheta(CScaling self, PetscReal *theta) {
  *theta = self->theta;
  return 0;
}
