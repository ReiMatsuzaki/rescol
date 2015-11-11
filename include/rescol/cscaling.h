#ifndef CSCALING_H
#define CSCALING_H
#ifdef __cplusplus
extern "C" {
#endif
#include <petscpf.h> 
/*
  arbitaly complex scaling
*/

typedef PF CScaling;
PetscErrorCode CScalingCreate(MPI_Comm comm, CScaling *p_self);
PetscErrorCode CScalingSetNone(CScaling self);
PetscErrorCode CScalingSetUniformCS(CScaling self, PetscReal t);
PetscErrorCode CScalingSetSharpECS(CScaling self, PetscReal r0, PetscReal t);
PetscErrorCode CScalingCalc(CScaling self, PetscReal *xs, int n,
			    PetscScalar *qrs, PetscScalar *Rrs);


#ifdef __cplusplus
}
#endif
#endif
