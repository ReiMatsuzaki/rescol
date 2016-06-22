#ifndef CSCALING_H
#define CSCALING_H
#ifdef __cplusplus
extern "C" {
#endif
#include <petscpf.h> 
/*
  arbitaly complex scaling
*/

struct _p_CScaling {
  MPI_Comm comm;
  PF pf;
  PetscBool use_cscaling;
  PetscReal R0;
  PetscReal theta;
} _p_CScaling;
typedef struct _p_CScaling* CScaling;

  // ---- Basics ----  
PetscErrorCode CScalingCreate(MPI_Comm comm, CScaling *p_self);
PetscErrorCode CScalingDestroy(CScaling *p_self);
PetscErrorCode CScalingView(CScaling slef, PetscViewer v);
PetscErrorCode CScalingSetNone(CScaling self);
PetscErrorCode CScalingSetUniformCS(CScaling self, PetscReal t);
PetscErrorCode CScalingSetSharpECS(CScaling self, PetscReal r0, PetscReal t);
PetscErrorCode CScalingSetFromOptions(CScaling self);

  // ---- Calculation ----
PetscErrorCode CScalingCalc(CScaling self, PetscReal *xs, int n,
			    PetscScalar *qrs, PetscScalar *Rrs);
PetscErrorCode CScalingCalcOne(CScaling self, PetscReal x,
			       PetscScalar *qr, PetscScalar *Rr);
PetscErrorCode CScalingQ(CScaling self, PetscBool *_use_cscaling);
PetscErrorCode CScalingGetRadius(CScaling self, PetscReal *R0);
PetscErrorCode CScalingGetTheta(CScaling self, PetscReal *theta);


#ifdef __cplusplus
}
#endif
#endif
