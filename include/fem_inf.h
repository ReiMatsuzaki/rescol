#ifndef FEM_INF_H
#define FEM_INF_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <petscmat.h>
#include "fd.h"
#include "bspline.h"
#include "dvr.h"
#include "viewerfunc.h"
#include "pot.h"

static const int FEMInfType_FD = 0;
static const int FEMInfType_BSS = 1;
static const int FEMInfType_DVR = 2;

// vtbl for FEM interface
typedef struct {
  PetscErrorCode (*Create)();
  PetscErrorCode (*Copy)();
  PetscErrorCode (*Destory)();
  PetscErrorCode (*View)();
  PetscErrorCode (*SetFromOptions)();

  PetscErrorCode (*Psi)();
  PetscErrorCode (*DerivPsi)();
  //PetscErrorCode (*GuessHEig)();
  PetscErrorCode (*GetSize)();

  PetscErrorCode (*SR1Mat)();
  PetscErrorCode (*D2R1Mat)();
  //  PetscErrorCode (*R2invR1Mat)();
  //  PetscErrorCode (*ENR1Mat)();
  PetscErrorCode (*PotR1Mat)();
  PetscErrorCode (*PotR1Vec)();
  PetscErrorCode (*EER2Mat)();
  PetscBool overlap_is_id;  
} FEMSc;

struct _p_FEMInf{
  MPI_Comm comm;
  FEMSc* sc; // scheme for polymorphism
  void* obj;  // address of object
  int type;
} ;
typedef struct _p_FEMInf* FEMInf;

// ----- getter of interface -----
PetscErrorCode FEMInfSetFD(FEMInf self, FD target);
PetscErrorCode FEMInfSetBSS(FEMInf self, BSS target);
PetscErrorCode FEMInfSetDVR(FEMInf self, DVR target);

// ---- Basic Method ------
PetscErrorCode FEMInfCreate(MPI_Comm comm, FEMInf *inf);
PetscErrorCode FEMInfDuplicate(FEMInf self, FEMInf *new_fem);  
PetscErrorCode FEMInfCopy(FEMInf self, FEMInf b);
PetscErrorCode FEMInfDestroy(FEMInf *inf);
//PetscErrorCode FEMInfCreateFD(FEMInf *inf, FD self);
PetscErrorCode FEMInfView(FEMInf self, PetscViewer v);
PetscErrorCode FEMInfViewFunc(FEMInf self, Vec c, ViewerFunc v);

// ---- Accessor ----
PetscErrorCode FEMInfSetFromOptions(FEMInf self);

PetscErrorCode FEMInfGetSize(FEMInf self, int *n);
PetscErrorCode FEMInfGetOverlapIsId(FEMInf self, PetscBool *is_id);

// ---- calculation -----
PetscErrorCode FEMInfPsi(FEMInf self, Vec cs, Vec xs, Vec ys);
PetscErrorCode FEMInfPsiOne(FEMInf self, Vec cs, PetscScalar x, PetscScalar *y);
PetscErrorCode FEMInfDerivPsi(FEMInf self, Vec c, Vec x, Vec y);
PetscErrorCode FEMInfDerivPsiOne(FEMInf self, Vec c, PetscReal x, PetscScalar *y);
PetscErrorCode FEMInfFit(FEMInf self, PF pf, KSP ksp, Vec c);
  //PetscErrorCode FEMInfGuessHEig(FEMInf self, int n, int l, PetscScalar z, Vec *v);

PetscErrorCode FEMInfCreateMat(FEMInf self, int dim, Mat *M);
PetscErrorCode FEMInfCreateVec(FEMInf self, int dim, Vec *v);

PetscErrorCode FEMInfSR1Mat(FEMInf self, Mat M);
  //PetscErrorCode FEMInfSR1MatNullable(FEMInf self, Mat M);
PetscErrorCode FEMInfD2R1Mat(FEMInf self, Mat M);
  //PetscErrorCode FEMInfR2invR1Mat(FEMInf self, Mat M);
  //PetscErrorCode FEMInfENR1Mat(FEMInf self, int q, double a, Mat M); 
PetscErrorCode FEMInfPotR1Mat(FEMInf self, Pot pot, Mat M);
PetscErrorCode FEMInfPotR1Vec(FEMInf self, Pot pot, Vec V);
PetscErrorCode FEMInfEER2Mat(FEMInf self, int q, Mat M); 

#ifdef __cplusplus
}
#endif
#endif
