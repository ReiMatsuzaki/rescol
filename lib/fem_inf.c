#include "../include/fem_inf.h"
#include <petscdraw.h>

// ----- Interface -----
FEMSc FD_Sc;
PetscErrorCode FEMInfSetFD(FEMInf self, FD target) {

  static int init = 0;
  if(init == 0) {
    FD_Sc.Create = FDCreate;
    FD_Sc.Copy = NULL;
    FD_Sc.Destory = FDDestroy;
    
    FD_Sc.View = FDView;
    FD_Sc.SetFromOptions = FDSetFromOptions;

    FD_Sc.Psi = NULL;
    FD_Sc.DerivPsi = NULL;
    //    FD_Sc.GuessHEig = FDGuessHEig;
    FD_Sc.GetSize = FDGetSize;
    
    FD_Sc.SR1Mat = FDSR1Mat;
    FD_Sc.D1R1Mat = NULL;
    FD_Sc.D2R1Mat = FDD2R1Mat;
    //    FD_Sc.R2invR1Mat = FDR2invR1Mat;
    //    FD_Sc.ENR1Mat = FDENR1Mat;    
    FD_Sc.PotR1Mat = NULL;
    FD_Sc.PotR1Vec = NULL;
    FD_Sc.EER2Mat = FDEER2Mat;        

    FD_Sc.overlap_is_id = PETSC_TRUE;    
    init = 1;
  }

  self->sc = &FD_Sc;
  self->obj = target;
  self->type = FEMInfType_FD;

  return 0;
}
FEMSc BSS_Sc;
PetscErrorCode FEMInfSetBSS(FEMInf self, BSS target) {

  static int init = 0;
  if(init == 0) {
    BSS_Sc.Create = BSSCreate;
    BSS_Sc.Copy = NULL;
    BSS_Sc.Destory = BSSDestroy;
    BSS_Sc.View = BSSView;
    BSS_Sc.SetFromOptions = BSSSetFromOptions;

    BSS_Sc.Psi = BSSPsi;
    BSS_Sc.DerivPsi = BSSDerivPsi;
    //    BSS_Sc.GuessHEig = NULL;
    BSS_Sc.GetSize = BSSGetSize;
    
    BSS_Sc.SR1Mat = BSSSR1Mat;
    BSS_Sc.D1R1Mat = NULL;
    BSS_Sc.D2R1Mat = BSSD2R1Mat;
    //    BSS_Sc.R2invR1Mat = BSSR2invR1Mat;
    //   BSS_Sc.ENR1Mat = BSSENR1Mat;    
    BSS_Sc.PotR1Mat = BSSPotR1Mat;    
    BSS_Sc.PotR1Vec = BSSPotR1Vec;
    BSS_Sc.EER2Mat = BSSEER2Mat;
    BSS_Sc.overlap_is_id = PETSC_FALSE;
    init = 1;
  }
  
  self->sc = &BSS_Sc;
  self->obj = target;
  self->type = FEMInfType_BSS;

  return 0;
}
FEMSc DVR_Sc;
PetscErrorCode FEMInfSetDVR(FEMInf self, DVR target) {

  static int init = 0;
  if(init == 0) {
    DVR_Sc.Create = DVRCreate;
    DVR_Sc.Copy = DVRCopy;
    DVR_Sc.Destory = DVRDestroy;
    
    DVR_Sc.View = DVRView;
    DVR_Sc.SetFromOptions = DVRSetFromOptions;

    DVR_Sc.Psi = DVRPsi;
    DVR_Sc.DerivPsi = DVRDerivPsi;
    DVR_Sc.GetSize = DVRGetSize;
    //    DVR_Sc.GuessHEig = NULL;

    DVR_Sc.SR1Mat = DVRSR1Mat;
    DVR_Sc.D1R1Mat = DVRD1R1Mat;
    DVR_Sc.D2R1Mat = DVRD2R1Mat;
    //    DVR_Sc.R2invR1Mat = DVRR2invR1Mat;
    //    DVR_Sc.ENR1Mat = DVRENR1Mat;
    DVR_Sc.PotR1Mat = DVRPotR1Mat;
    DVR_Sc.PotR1Vec = DVRPotR1Vec;
    DVR_Sc.EER2Mat = DVREER2Mat;
    DVR_Sc.overlap_is_id = PETSC_TRUE;
    init = 1;
  }
  
  self->sc = &DVR_Sc;
  self->obj = target;
  self->type = FEMInfType_DVR;
  
  return 0;
}

// ---- Basic Method ------
PetscErrorCode FEMInfCreate(MPI_Comm comm, FEMInf *p_self) {
  FEMInf self;
  PetscNew(&self);
  
  self->comm = comm;
  self->sc = NULL;
  self->obj = NULL;
  self->type = -1;

  *p_self = self;
  return 0;
}
PetscErrorCode FEMInfDuplicate(FEMInf self, FEMInf *new_fem) {
  //  if(self->sc->Duplicate == NULL)
  //    SETERRQ(self->comm, 1, "method is null");

  FEMInfCreate(self->comm, new_fem);
  if(self->type == FEMInfType_FD) {
    FD ptr;
    FDCreate(self->comm, &ptr);
    FEMInfSetFD(*new_fem, ptr);
  }
  if(self->type == FEMInfType_BSS) {
    BSS ptr;
    BSSCreate(self->comm, &ptr);
    FEMInfSetBSS(*new_fem, ptr);
  }
  if(self->type == FEMInfType_DVR) {
    DVR ptr;
    DVRCreate(self->comm, &ptr);
    FEMInfSetDVR(*new_fem, ptr);
  }
  return 0;
}
PetscErrorCode FEMInfCopy(FEMInf self, FEMInf other) {

  if(self->sc->Copy == NULL)
    SETERRQ(self->comm, 1, "method is null");
  if(self->type != other->type) 
    SETERRQ(self->comm, 1, "self and other is not same object");
  
  PetscErrorCode ierr;
  ierr = self->sc->Copy(self->obj, other->obj); CHKERRQ(ierr);
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

  PetscViewerASCIIPrintf(v, "FEMInf object:\n");  
  PetscViewerASCIIPushTab(v);
  self->sc->View(self->obj, v);
  PetscViewerASCIIPopTab(v);

  return 0;
}
PetscErrorCode FEMInfViewFunc_ASCII(FEMInf self, Vec c, ViewerFunc v) {
    PetscInt num;
    PetscReal *xs;
    ViewerFuncGetRange(v, &num, &xs);

#if defined(PETSC_USE_COMPLEX)
    PetscViewerASCIIPrintf(v->base, "r,re_y,im_y\n");
#else
      PetscViewerASCIIPrintf(v->base, "r,y\n"); CHKERRQ(ierr);
#endif    

    for(int i = 0; i < num; i++) {
      PetscReal x = xs[i];
      PetscScalar y; FEMInfPsiOne(self, c, x, &y);
#if defined(PETSC_USE_COMPLEX)
      PetscReal re = PetscRealPart(y);
      PetscReal im = PetscImaginaryPart(y);
      PetscViewerASCIIPrintf(v->base, "%f,%f,%f\n", x, re, im);
#else
      PetscViewerASCIIPrintf(v->base, "%f,%f\n", x, y); CHKERRQ(ierr);
#endif
    }
    
    return 0;
    
}
PetscErrorCode FEMInfViewFunc_Draw(FEMInf self, Vec c, ViewerFunc v) {

  PetscDraw draw; PetscViewerDrawGetDraw(v->base, 0, &draw);
  PetscDrawLG lg; PetscDrawLGCreate(draw, 1, &lg);

  PetscInt num;
  PetscReal *xs;
  ViewerFuncGetRange(v, &num, &xs);

  for(int i = 0; i < num; i++) {
    PetscReal x = xs[i];
    PetscScalar y[1]; FEMInfPsiOne(self, c, x, &y[0]);
#if defined(PETSC_USE_COMPLEX)
    PetscReal re[1] = {PetscRealPart(y[0])};
    //    PetscReal im[1] = {PetscImaginaryPart(y[0])};
    PetscDrawLGAddPoint(lg, &xs[i], re);
#else
    PetscDrawLGAddPoint(lg, &xs[i], y);
#endif
  }  
  PetscDrawLGDraw(lg);
  return 0;
}
PetscErrorCode FEMInfViewFunc(FEMInf self, Vec c, ViewerFunc v) {

  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  //  PetscViewerType type;     PetscViewerGetType(v, &type);
  //  PetscViewerFormat format; PetscViewerGetFormat(v, &format);

  ierr = PetscObjectTypeCompare((PetscObject)v->base,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v->base,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v->base,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);    
  if(iascii) {
    ierr = FEMInfViewFunc_ASCII(self, c, v); CHKERRQ(ierr);
  } else if(isdraw) {
    ierr = FEMInfViewFunc_Draw(self, c, v); CHKERRQ(ierr);
  } else {
  }
  return 0;
}

// ---- Accessor ----
PetscErrorCode FEMInfSetFromOptions(FEMInf self) {

  char type[10];
  PetscBool find;
  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, NULL, "-fem_type", type, 10, &find); CHKERRQ(ierr);
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

// ---- Calculation ----
PetscErrorCode FEMInfFit(FEMInf self, PF pf, KSP ksp, Vec c) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfFit", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
  PetscErrorCode ierr;
  PetscBool is_id; FEMInfGetOverlapIsId(self, &is_id);

  if(is_id) {
    ierr = FEMInfPotR1Vec(self, pf, c); CHKERRQ(ierr);
  } else {
    if(ksp == NULL) {
      SETERRQ(self->comm, 1, "ksp is null");
    }
    Vec V; FEMInfCreateVec(self, 1, &V); FEMInfPotR1Vec(self, pf, V);
    Mat S; FEMInfCreateMat(self, 1, &S); FEMInfSR1Mat(self, S);
    ierr = KSPSetOperators(ksp, S, S); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, V, c); CHKERRQ(ierr);
    MatDestroy(&S);
    VecDestroy(&V);
  }

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
    
  return 0;
  
}
PetscErrorCode FEMInfGetOverlapIsId(FEMInf this, PetscBool *is_id) {
  *is_id = this->sc->overlap_is_id;  
  return 0;
}
PetscErrorCode FEMInfGetSize(FEMInf this, int *n) {
  
  if(this->sc->GetSize == NULL)
    SETERRQ(this->comm, 1, "method is null: GetSize");

  this->sc->GetSize(this->obj, n);
  return 0;
}

// ---- Calculation ----
PetscErrorCode FEMInfPsi(FEMInf self, Vec cs, Vec xs, Vec ys) {
  /*
    gives function values on x. { f(x) | x in xs}
    cs: coefficient vector
    xs: x list
    ys: result
   */

  int EVENT_id;
  PetscLogEventRegister("FEMInfPsi", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  
  
  if(self->sc->Psi == NULL) {
    SETERRQ(self->comm, 1, "method is null: Psi");
  }

  PetscErrorCode ierr;
  ierr = self->sc->Psi(self->obj, cs, xs, ys); CHKERRQ(ierr);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;
}
PetscErrorCode FEMInfPsiOne(FEMInf self, Vec cs, PetscScalar x, PetscScalar *y) {
  /*
    gives function values on x. { f(x) | x in xs}
    cs: coefficient vector
    x: x
    y: result
   */
  PetscErrorCode ierr;
  Vec xs;
  ierr = VecCreate(self->comm, &xs); CHKERRQ(ierr);
  ierr = VecSetSizes(xs, PETSC_DECIDE, 1); CHKERRQ(ierr);
  ierr = VecSetUp(xs);
  ierr = VecSetValue(xs, 0, x, INSERT_VALUES); CHKERRQ(ierr);
  Vec ys;
  ierr = VecDuplicate(xs, &ys); CHKERRQ(ierr);

  ierr = self->sc->Psi(self->obj, cs, xs, ys); CHKERRQ(ierr);

  PetscScalar *y_ptr;
  ierr = VecGetArray(ys, &y_ptr); CHKERRQ(ierr);
  *y = y_ptr[0];
  ierr = VecRestoreArray(ys, &y_ptr); CHKERRQ(ierr);

  VecDestroy(&xs);
  VecDestroy(&ys);

  return 0;
}
PetscErrorCode FEMInfDerivPsi(FEMInf self, Vec c, Vec xs, Vec ys) {

  if(self->sc->DerivPsi == NULL)
    SETERRQ(self->comm, 1, "method is null: DerivPsi");
  
  self->sc->DerivPsi(self->obj, c, xs, ys);
  return 0;

}
PetscErrorCode FEMInfDerivPsiOne(FEMInf self, Vec c, PetscReal x, PetscScalar *y) {

  PetscErrorCode ierr;
  Vec xs;
  ierr = VecCreate(self->comm, &xs); CHKERRQ(ierr);
  ierr = VecSetSizes(xs, PETSC_DECIDE, 1); CHKERRQ(ierr);
  ierr = VecSetUp(xs);
  ierr = VecSetValue(xs, 0, x, INSERT_VALUES); CHKERRQ(ierr);
  Vec ys;
  ierr = VecDuplicate(xs, &ys); CHKERRQ(ierr);

  ierr = FEMInfDerivPsi(self, c, xs, ys);

  PetscScalar *y_ptr;
  ierr = VecGetArray(ys, &y_ptr); CHKERRQ(ierr);
  *y = y_ptr[0];
  ierr = VecRestoreArray(ys, &y_ptr); CHKERRQ(ierr);

  VecDestroy(&xs);
  VecDestroy(&ys);

  return 0;  

}
// - to be removed
PetscErrorCode FEMInfGuessHEig(FEMInf self, int n, int l, PetscScalar z, Vec *v) {


  SETERRQ(self->comm, 1, "unsupported method");
  /*
  if(self->sc->GuessHEig == NULL)
    SETERRQ(self->comm, 1, "method is null: GuessHEig");

  self->sc->GuessHEig(self->obj, n, l, z, v);
  return 0;
  */
}

PetscErrorCode FEMInfCreateMat(FEMInf self, int dim, Mat *M) {
  if(dim < 1)
    SETERRQ(self->comm, 1, "dim must be positive integer");

  int n; FEMInfGetSize(self, &n);
  int nn = pow(n, dim);

  if(n <= 0) {
    printf("size of basis = %d\n", n);
    SETERRQ(self->comm, 1, "size of basis is negative or 0.");
  }

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
    SETERRQ(self->comm, 1, "method is null: SR1Mat");

  int EVENT_id;
  PetscLogEventRegister("FEMInfSR1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  
  self->sc->SR1Mat(self->obj, M);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;
}
// - to be removed
PetscErrorCode FEMInfSR1MatNullable(FEMInf self, Mat M) {
  SETERRQ(self->comm, 1, "unsupported method");
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
PetscErrorCode FEMInfD1R1Mat(FEMInf self, Mat M) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfD1R1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  
  if(self->sc->D1R1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: D1R1Mat");

  self->sc->D1R1Mat(self->obj, M);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;

}
PetscErrorCode FEMInfD2R1Mat(FEMInf self, Mat M) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfD2R1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
  if(self->sc->D2R1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: D2R1Mat");

  self->sc->D2R1Mat(self->obj, M);


  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;

}
// - to be removed
PetscErrorCode FEMInfR2invR1Mat(FEMInf self, Mat M) {
  SETERRQ(self->comm, 1, "unsupported method");
  /*
  if(self->sc->R2invR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: R2invR1Mat");

  self->sc->R2invR1Mat(self->obj, M);
  return 0;
  */
}
// - to be removed
PetscErrorCode FEMInfENR1Mat(FEMInf self, int q, double a, Mat M) {
  SETERRQ(self->comm, 1, "unsupported method");
  /*
  if(self->sc->ENR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: ENR1Mat");

  self->sc->ENR1Mat(self->obj, q, a, M);
  */
  return 0;

}
PetscErrorCode FEMInfPotR1Mat(FEMInf self, Pot pot, Mat M) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfPotR1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
  if(self->sc->PotR1Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: PotR1Mat");

  self->sc->PotR1Mat(self->obj, pot, M);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;
}
PetscErrorCode FEMInfPotR1Vec(FEMInf self, Pot pot, Vec V) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfPotR1Vec", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
  if(self->sc->PotR1Vec == NULL)
    SETERRQ(self->comm, 1, "method is null: PotR1Vec");

  self->sc->PotR1Vec(self->obj, pot, V);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);

  return 0;

}
PetscErrorCode FEMInfEER2Mat(FEMInf self, int q, Mat M) {

  int EVENT_id;
  PetscLogEventRegister("FEMInfEER2Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
  if(self->sc->EER2Mat == NULL)
    SETERRQ(self->comm, 1, "method is null: EER2Mat");

  self->sc->EER2Mat(self->obj, q, M);


  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;

}

