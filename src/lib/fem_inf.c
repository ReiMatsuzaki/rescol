#include <rescol/fem_inf.h>
#include <rescol/dvr.h>
#include <rescol/bspline.h>

// ----- Interface -----
PetscErrorCode FEMInfCreate(FEMInf *inf, MPI_Comm comm) {
  FEMInf _inf;
  PetscNew(&_inf);
  
  _inf->comm = comm;
  _inf->sc = NULL;
  _inf->obj = NULL;

  *inf = _inf;
  return 0;
}
FEMSc FD_Sc;
PetscErrorCode FEMInfCreateFD(FEMInf *inf, FD this) {

  FEMInfCreate(inf, this->comm);

  static int init = 0;
  if(init == 0) {
    FD_Sc.Create = FDCreate;
    FD_Sc.CreateFromOptions = FDCreateFromOptions;
    FD_Sc.Destory = FDDestroy;
    FD_Sc.FPrintf = FDFPrintf;
    FD_Sc.SetSR1Mat = FDSetSR1Mat;
    FD_Sc.SetD2R1Mat = FDSetD2R1Mat;
    FD_Sc.SetR2invR1Mat = FDSetR2invR1Mat;
    FD_Sc.SetENR1Mat = FDSetENR1Mat;    
    FD_Sc.SetEER2Mat = FDSetEER2Mat;        
    FD_Sc.BasisPsi = NULL;
    FD_Sc.GuessHEig = FDGuessHEig;
    FD_Sc.GetSize = FDGetSize;
    FD_Sc.overlap_is_id = PETSC_TRUE;    
    init = 1;
  }

  (*inf)->sc = &FD_Sc;
  (*inf)->obj = this;

  return 0;
}
FEMSc BSS_Sc;
PetscErrorCode FEMInfCreateBSS(FEMInf *inf, BSS this) {

  FEMInfCreate(inf, this->comm);

  static int init = 0;
  if(init == 0) {
    BSS_Sc.Create = BSSCreate;
    BSS_Sc.CreateFromOptions = BSSCreateFromOptions;
    BSS_Sc.Destory = BSSDestroy;
    BSS_Sc.FPrintf = BSSFPrintf;
    BSS_Sc.SetSR1Mat = BSSSetSR1Mat;
    BSS_Sc.SetD2R1Mat = BSSSetD2R1Mat;
    BSS_Sc.SetR2invR1Mat = BSSSetR2invR1Mat;
    BSS_Sc.SetENR1Mat = BSSSetENR1Mat;    
    BSS_Sc.SetEER2Mat = BSSSetEER2Mat;
    BSS_Sc.BasisPsi = BSSBasisPsi;
    BSS_Sc.GuessHEig = NULL;
    BSS_Sc.GetSize = BSSGetSize;
    BSS_Sc.overlap_is_id = PETSC_FALSE;
    init = 1;
  }
  
  (*inf)->sc = &BSS_Sc;
  (*inf)->obj = this;
  
  return 0;
}
FEMSc DVR_Sc;
PetscErrorCode FEMInfCreateDVR(FEMInf *inf, DVR this) {

  FEMInfCreate(inf, this->comm);

  static int init = 0;
  if(init == 0) {
    DVR_Sc.Create = DVRCreate;
    DVR_Sc.CreateFromOptions = DVRCreateFromOptions;
    DVR_Sc.Destory = DVRDestroy;
    DVR_Sc.FPrintf = DVRFPrintf;
    DVR_Sc.SetSR1Mat = DVRSetSR1Mat;
    DVR_Sc.SetD2R1Mat = DVRSetD2R1Mat;
    DVR_Sc.SetR2invR1Mat = DVRSetR2invR1Mat;
    DVR_Sc.SetENR1Mat = DVRSetENR1Mat;    
    DVR_Sc.SetEER2Mat = DVRSetEER2Mat;        
    DVR_Sc.BasisPsi = DVRBasisPsi;
    DVR_Sc.GetSize = DVRGetSize;
    DVR_Sc.GuessHEig = NULL;
    DVR_Sc.overlap_is_id = PETSC_TRUE;
    init = 1;
  }
  
  (*inf)->sc = &DVR_Sc;
  (*inf)->obj = this;
  
  return 0;
}

// ---- Basic Method ------
PetscErrorCode FEMInfCreateFromOptions(FEMInf *inf, MPI_Comm comm) {

  char type[10];
  PetscBool find;
  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, "-fem_type", type, 10, &find); CHKERRQ(ierr);
  if (strcmp(type, "fd") == 0){
    FD fd;
    ierr = FDCreateFromOptions(&fd, comm);
    ierr = FEMInfCreateFD(inf, fd);

  }else if(strcmp(type, "bss") == 0) {
    BSS bss;
    ierr = BSSCreateFromOptions(&bss, comm);
    ierr = FEMInfCreateBSS(inf, bss);
  } else if(strcmp(type, "dvr") == 0) {
    DVR dvr;
    ierr = DVRCreateFromOptions(&dvr, comm); CHKERRQ(ierr);
    ierr = FEMInfCreateDVR(inf, dvr);
  } else {
    SETERRQ(comm, 1, "-fem_type <- {bss, dvr}");
  }
  
  return 0;
}
PetscErrorCode FEMInfDestroy(FEMInf *inf) {

  if((*inf)->sc->Destory == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");
  
  (*inf)->sc->Destory(&((*inf)->obj));

  PetscFree(*inf);

  return 0;
}
PetscErrorCode FEMInfFPrintf(FEMInf this, FILE *file, int lvl) {

  if(this->sc->FPrintf == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->FPrintf(this->obj, file, lvl);
  return 0;

}
PetscErrorCode FEMInfView(FEMInf this) {
  FEMInfFPrintf(this, stdout, 0);
  return 0;
}

// ---- Accessor ----
PetscErrorCode FEMInfGetOverlapIsId(FEMInf this, PetscBool *is_id) {
  *is_id = this->sc->overlap_is_id;  
  return 0;
}
PetscErrorCode FEMInfGetSize(FEMInf this, int *n) {
  
  if(this->sc->GetSize == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->GetSize(this->obj, n);
  return 0;
}

// ---- Calculation ----
PetscErrorCode FEMInfSetSR1Mat(FEMInf this, Mat *M) {

  if(this->sc->SetSR1Mat == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->SetSR1Mat(this->obj, M);
  return 0;
}
PetscErrorCode FEMInfSetD2R1Mat(FEMInf this, Mat *M) {

  if(this->sc->SetD2R1Mat == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->SetD2R1Mat(this->obj, M);
  return 0;

}
PetscErrorCode FEMInfSetR2invR1Mat(FEMInf this, Mat *M) {

  if(this->sc->SetR2invR1Mat == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->SetR2invR1Mat(this->obj, M);
  return 0;
}
PetscErrorCode FEMInfSetENR1Mat(FEMInf this, int q, double a, Mat *M) {

  if(this->sc->SetENR1Mat == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->SetENR1Mat(this->obj, q, a, M);
  return 0;

}
PetscErrorCode FEMInfSetEER2Mat(FEMInf this, int q, Mat *M) {

  if(this->sc->SetEER2Mat == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->SetEER2Mat(this->obj, q, M);
  return 0;

}
PetscErrorCode FEMInfBasisPsi(FEMInf this, int i, PetscScalar x, PetscScalar *y) {

  if(this->sc->BasisPsi == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->BasisPsi(this->obj, i, x, y);
  return 0;
}

PetscErrorCode FEMInfGuessHEig(FEMInf this, int n, int l, PetscScalar z, Vec *v) {

  if(this->sc->GuessHEig == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "method is null");

  this->sc->GuessHEig(this->obj, n, l, z, v);
  return 0;
}

