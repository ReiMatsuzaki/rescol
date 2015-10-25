#ifndef FEM_INF_H
#define FEM_INF_H

#include <petscmat.h>
#include "fd.h"
#include "bspline.h"
#include "dvr.h"


typedef struct {
  PetscErrorCode (*Create)();
  PetscErrorCode (*CreateFromOptions)();
  PetscErrorCode (*Destory)();
  PetscErrorCode (*FPrintf)();
  PetscErrorCode (*SetSR1Mat)();
  PetscErrorCode (*SetD2R1Mat)();
  PetscErrorCode (*SetR2invR1Mat)();
  PetscErrorCode (*SetENR1Mat)();
  PetscErrorCode (*SetEER2Mat)();
  PetscErrorCode (*BasisPsi)();
  PetscErrorCode (*GuessHEig)();
  PetscBool overlap_is_id;
  
} FEMSc;

typedef struct {
  FEMSc* sc; // scheme for polymorphism
  void* obj;  // address of object
} FEMInf;

// ----- getter of interface -----
PetscErrorCode FEMInfCreateFD(FEMInf *inf, FD this);
PetscErrorCode FEMInfCreateBSS(FEMInf *inf, BSS this);
PetscErrorCode FEMInfCreateDVR(FEMInf *inf, DVR this);

// ---- method ------
PetscErrorCode FEMInfCreateFromOptions(FEMInf *this, MPI_Comm comm);
PetscErrorCode FEMInfDestroy(FEMInf *inf);
PetscErrorCode FEMInfFPrintf(FEMInf this, FILE *file, int lvl);
PetscErrorCode FEMInfSetSR1Mat(FEMInf this, Mat *M);
PetscErrorCode FEMInfSetD2R1Mat(FEMInf this, Mat *M);
PetscErrorCode FEMInfSetR2invR1Mat(FEMInf this, Mat *M);
PetscErrorCode FEMInfSetENR1Mat(FEMInf this, int q, double a, Mat *M); 
PetscErrorCode FEMInfSetEER2Mat(FEMInf this, int q, Mat *M); 
PetscErrorCode FEMInfBasisPsi(FEMInf this, int i, PetscScalar x, PetscScalar *y);
PetscErrorCode FEMInfGetOverlapIsId(FEMInf this, PetscBool *is_id);
PetscErrorCode FEMInfGuessHEig(FEMInf this, int n, int l, Vec *v);

#endif
