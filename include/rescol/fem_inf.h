#ifndef FEM_INF_H
#define FEM_INF_H

#ifdef __cplusplus
extern "C" {
#endif 
#include <petscmat.h>
#include "fd.h"
#include "bspline.h"
#include "dvr.h"
#include "pot.h"

// vtbl for FEM interface
typedef struct {
  PetscErrorCode (*Create)();
  PetscErrorCode (*CreateFromOptions)();
  PetscErrorCode (*Destory)();
  PetscErrorCode (*FPrintf)();
  PetscErrorCode (*SetSR1Mat)();
  PetscErrorCode (*SetD2R1Mat)();
  PetscErrorCode (*SetR2invR1Mat)();
  PetscErrorCode (*SetENR1Mat)();
  PetscErrorCode (*SetPotR1Mat)();
  PetscErrorCode (*SetEER2Mat)();
  PetscErrorCode (*BasisPsi)();
  //  PetscErrorCode (*Psi)();
  PetscErrorCode (*GuessHEig)();
  PetscErrorCode (*GetSize)();
  PetscBool overlap_is_id;
  
} FEMSc;

struct _p_FEMInf{
  MPI_Comm comm;
  FEMSc* sc; // scheme for polymorphism
  void* obj;  // address of object
} ;
typedef struct _p_FEMInf* FEMInf;

// ----- getter of interface -----
PetscErrorCode FEMInfCreateFD(FEMInf *inf, FD self);
PetscErrorCode FEMInfCreateBSS(FEMInf *inf, BSS self);
PetscErrorCode FEMInfCreateDVR(FEMInf *inf, DVR self);

// ---- Basic Method ------
PetscErrorCode FEMInfCreateFromOptions(FEMInf *self, MPI_Comm comm);
PetscErrorCode FEMInfDestroy(FEMInf *inf);
PetscErrorCode FEMInfFPrintf(FEMInf self, FILE *file, int lvl);
PetscErrorCode FEMInfView(FEMInf self);

// ----- Accessor ----
PetscErrorCode FEMInfGetSize(FEMInf self, int *n);
PetscErrorCode FEMInfGetOverlapIsId(FEMInf self, PetscBool *is_id);

// ---- calculation -----
PetscErrorCode FEMInfSetSR1Mat(FEMInf self, Mat *M);
PetscErrorCode FEMInfSEtSR1MatNullable(FEMInf self, Mat *M);
PetscErrorCode FEMInfSetD2R1Mat(FEMInf self, Mat *M);
PetscErrorCode FEMInfSetR2invR1Mat(FEMInf self, Mat *M);
PetscErrorCode FEMInfSetENR1Mat(FEMInf self, int q, double a, Mat *M); 
PetscErrorCode FEMInfSetPOTR1Mat(FEMInf self, POT pot, Mat *M);
PetscErrorCode FEMInfSetEER2Mat(FEMInf self, int q, Mat *M); 
PetscErrorCode FEMInfBasisPsi(FEMInf self, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode FEMInfPsi(FEMInf self, PetscReal x, Vec c, PetscScalar *y);
/*
PetscErrorCode FEMInfWritePsi(FEMInf self, PetscReal *xs, int num_x, 
			      Vec cs, FILE *file);
*/
PetscErrorCode FEMInfGuessHEig(FEMInf self, int n, int l, PetscScalar z, Vec *v);

#ifdef __cplusplus
}
#endif
#endif
