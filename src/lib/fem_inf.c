#include <rescol/fem_inf.h>

// ----- Interface -----

FEMSc FD_Sc;
PetscErrorCode FEMInfSetFD(FEMInf self, FD target) {

  static int init = 0;
  if(init == 0) {
    FD_Sc.Create = FDCreate;
    FD_Sc.Destory = FDDestroy;
    FD_Sc.View = FDView;
    FD_Sc.SetFromOptions = FDSetFromOptions;

    FD_Sc.Psi = NULL;
    FD_Sc.GuessHEig = FDGuessHEig;
    FD_Sc.GetSize = FDGetSize;
    
    FD_Sc.SR1Mat = FDSR1Mat;
    FD_Sc.D2R1Mat = FDD2R1Mat;
    FD_Sc.R2invR1Mat = FDR2invR1Mat;
    FD_Sc.ENR1Mat = FDENR1Mat;    
    FD_Sc.PotR1Mat = NULL;
    FD_Sc.EER2Mat = FDEER2Mat;        

    FD_Sc.overlap_is_id = PETSC_TRUE;    
    init = 1;
  }

  self->sc = &FD_Sc;
  self->obj = target;

  return 0;
}
FEMSc BSS_Sc;
PetscErrorCode FEMInfSetBSS(FEMInf self, BSS target) {

  static int init = 0;
  if(init == 0) {
    BSS_Sc.Create = BSSCreate;
    BSS_Sc.Destory = BSSDestroy;
    BSS_Sc.View = BSSView;
    BSS_Sc.SetFromOptions = BSSSetFromOptions;

    BSS_Sc.Psi = BSSPsi;
    BSS_Sc.GuessHEig = NULL;
    BSS_Sc.GetSize = BSSGetSize;
    
    BSS_Sc.SR1Mat = BSSSR1Mat;
    BSS_Sc.D2R1Mat = BSSD2R1Mat;
    BSS_Sc.R2invR1Mat = BSSR2invR1Mat;
    BSS_Sc.ENR1Mat = BSSENR1Mat;    
    BSS_Sc.PotR1Mat = BSSPotR1Mat;    
    BSS_Sc.EER2Mat = BSSEER2Mat;
    BSS_Sc.overlap_is_id = PETSC_FALSE;
    init = 1;
  }
  
  self->sc = &BSS_Sc;
  self->obj = target;
  
  return 0;
}
FEMSc DVR_Sc;
PetscErrorCode FEMInfSetDVR(FEMInf self, DVR target) {

  static int init = 0;
  if(init == 0) {
    DVR_Sc.Create = DVRCreate;
    DVR_Sc.Destory = DVRDestroy;
    DVR_Sc.View = DVRView;
    DVR_Sc.SetFromOptions = DVRSetFromOptions;

    DVR_Sc.Psi = NULL;
    DVR_Sc.GetSize = DVRGetSize;
    DVR_Sc.GuessHEig = NULL;

    DVR_Sc.SR1Mat = DVRSR1Mat;
    DVR_Sc.D2R1Mat = DVRD2R1Mat;
    DVR_Sc.R2invR1Mat = DVRR2invR1Mat;
    DVR_Sc.ENR1Mat = DVRENR1Mat;
    DVR_Sc.PotR1Mat = NULL;
    DVR_Sc.EER2Mat = NULL;
    DVR_Sc.overlap_is_id = PETSC_TRUE;
    init = 1;
  }
  
  self->sc = &DVR_Sc;
  self->obj = target;
  
  return 0;
}

// ---- Basic Method ------
PetscErrorCode FEMInfCreate(MPI_Comm comm, FEMInf *p_self) {
  FEMInf self;
  PetscNew(&self);
  
  self->comm = comm;
  self->sc = NULL;
  self->obj = NULL;

  *p_self = self;
  return 0;
}
PetscErrorCode FEMInfDestroy(FEMInf *inf) {

  PetscErrorCode ierr;

  if((*inf)->sc->Destory == NULL)
    SETERRQ((*inf)->comm, 1, "method is null");
  
  ierr = (*inf)->sc->Destory(&((*inf)->obj)); CHKERRQ(ierr);
  
  ierr = PetscFree(*inf); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode FEMInfView(FEMInf self, PetscViewer v) {
  if(self->sc->View == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->View(self->obj, v);

  return 0;
}
PetscErrorCode FEMInfViewFunc(FEMInf self, Vec c, ViewerFunc v) {

  PetscViewerType type;
  PetscViewerGetType(v->base, &type);
  if(strcmp(type, "ascii") != 0) {
    char msg[100]; sprintf(msg, "unsupported type: %s", type);
    SETERRQ(self->comm, 1, msg);
  }

  PetscInt num;
  PetscReal *xs;
  ViewerFuncGetXs(v, &num, &xs);

  for(int i = 0; i < num; i++) {
    PetscReal x = xs[i];
    PetscScalar y; FEMInfPsi(self, c, x, &y);
#if defined(PETSC_USE_COMPLEX)
    PetscReal re = PetscRealPart(y);
    PetscReal im = PetscImaginaryPart(y);
    PetscViewerASCIIPrintf(v->base, "%f %f %f\n", x, re, im);
#else
    PetscViewerASCIIPrintf(v->base, "%f %f\n", x, y); CHKERRQ(ierr);
#endif
  }
  return 0;
}

// ---- Accessor ----
PetscErrorCode FEMInfSetFromOptions(FEMInf self) {

  char type[10];
  PetscBool find;
  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, "-fem_type", type, 10, &find); CHKERRQ(ierr);
  if (strcmp(type, "fd") == 0){
    FD fd; FDCreate(self->comm, &fd);
    ierr = FDSetFromOptions(fd); CHKERRQ(ierr);
    ierr = FEMInfSetFD(self, fd); CHKERRQ(ierr);

  }else if(strcmp(type, "bss") == 0) {
    BSS bss; BSSCreate(self->comm, &bss);
    ierr = BSSSetFromOptions(bss); CHKERRQ(ierr);
    ierr = FEMInfSetBSS(self, bss); CHKERRQ(ierr);

  } else if(strcmp(type, "dvr") == 0) {
    DVR dvr; DVRCreate(self->comm, &dvr);
    ierr = DVRSetFromOptions(dvr); CHKERRQ(ierr);
    ierr = FEMInfSetDVR(self, dvr);
  } else if(strcmp(type, "") == 0){
    SETERRQ(self->comm, 1, 
	    "value of option -fem_type is empty. chose {fd, bss, dvr}.");

  } else {
    SETERRQ(self->comm, 1, "-fem_type <- {fd, bss, dvr}");
  }
  
  return 0;
}

PetscErrorCode FEMInfGetOverlapIsId(FEMInf this, PetscBool *is_id) {
  *is_id = this->sc->overlap_is_id;  
  return 0;
}
PetscErrorCode FEMInfGetSize(FEMInf this, int *n) {
  
  if(this->sc->GetSize == NULL)
    SETERRQ(this->comm, 1, "method is null");

  this->sc->GetSize(this->obj, n);
  return 0;
}

// ---- Calculation ----
PetscErrorCode FEMInfPsi(FEMInf self, Vec c, PetscReal x, PetscScalar *y) {

  if(self->sc->Psi == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->Psi(self->obj, c, x, y);
  return 0;
}
PetscErrorCode FEMInfGuessHEig(FEMInf self, int n, int l, PetscScalar z, Vec *v) {

  if(self->sc->GuessHEig == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->GuessHEig(self->obj, n, l, z, v);
  return 0;
}

PetscErrorCode FEMInfCreateMat(FEMInf self, int dim, Mat *M) {
  if(dim < 1)
    SETERRQ(self->comm, 1, "dim must be positive integer");

  int n; FEMInfGetSize(self, &n);
  int nn = pow(n, dim);

  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nn, nn);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode FEMInfCreateVec(FEMInf self, int dim, Vec *v) {
  if(dim < 1)
    SETERRQ(self->comm, 1, "dim must be positive integer");

  int n; FEMInfGetSize(self, &n);
  int nn = pow(n, dim);

  VecCreate(self->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, nn);
  VecSetUp(*v);
  return 0;
}

PetscErrorCode FEMInfSR1Mat(FEMInf self, Mat M) {

  if(self->sc->SR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->SR1Mat(self->obj, M);
  return 0;
}
PetscErrorCode FEMInfSR1MatNullable(FEMInf self, Mat M) {
  PetscErrorCode ierr;
  PetscBool s_is_id;
  ierr = FEMInfGetOverlapIsId(self, &s_is_id); CHKERRQ(ierr);
  if(s_is_id)
    M = NULL;
  else {
    ierr =FEMInfSR1Mat(self, M); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode FEMInfD2R1Mat(FEMInf self, Mat M) {

  if(self->sc->D2R1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->D2R1Mat(self->obj, M);
  return 0;

}
PetscErrorCode FEMInfR2invR1Mat(FEMInf self, Mat M) {

  if(self->sc->R2invR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->R2invR1Mat(self->obj, M);
  return 0;
}
PetscErrorCode FEMInfENR1Mat(FEMInf self, int q, double a, Mat M) {

  if(self->sc->ENR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->ENR1Mat(self->obj, q, a, M);
  return 0;

}
PetscErrorCode FEMInfPotR1Mat(FEMInf self, Pot pot, Mat M) {

  if(self->sc->PotR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->PotR1Mat(self->obj, pot, M);
  return 0;
}
PetscErrorCode FEMInfEER2Mat(FEMInf self, int q, Mat M) {

  if(self->sc->EER2Mat == NULL)
    SETERRQ(self->comm, 1, "method is null");

  self->sc->EER2Mat(self->obj, q, M);
  return 0;

}

