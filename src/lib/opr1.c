#include <rescol/opr1.h>

PetscErrorCode D2R1MatBSS(void *obj, BSS bss, Mat M) {
  PetscErrorCode ierr;
  ierr = BSSD2R1Mat(bss, M);
  return 0;
}
PetscErrorCode D2R1MatDVR(void *obj, DVR a, Mat M) {
  PetscErrorCode ierr;
  ierr = DVRD2R1Mat(a, M);
  return 0;
}
PetscErrorCode D2R1MatFD(void *obj, FD a, Mat M) {
  PetscErrorCode ierr;
  ierr = FDD2R1Mat(a, M);
  return 0;
}
PetscErrorCode D2R1MatView(void *obj, PetscViewer v) {

  PetscViewerASCIIPrintf(v, "name: D2\n");
  return 0;
}

PetscErrorCode PfR1MatBSS(void *obj, BSS bss, Mat M) {
  PetscErrorCode ierr;
  PF pot = (PF)obj;
  ierr = BSSPotR1Mat(bss, pot, M);
  return 0;
}
PetscErrorCode PfR1MatDVR(void *obj, DVR dvr, Mat M) {
  PetscErrorCode ierr;
  PF pot = (PF)obj;
  MPI_Comm comm; PetscObjectGetComm((PetscObject)pot, &comm);
  SETERRQ(comm, 1, "not implemented yet");
  //ierr = DVRPotR1Mat(bss, *pot, M);
  return 0;
}
PetscErrorCode PfR1MatFD(void *obj, FD dvr, Mat M) {
  PetscErrorCode ierr;
  PF pot = (PF)obj;
  MPI_Comm comm; PetscObjectGetComm((PetscObject)pot, &comm);
  SETERRQ(comm, 1, "not implemented yet");
  //ierr = FDPotR1Mat(bss, *pot, M);
  return 0;
}
PetscErrorCode PfR1MatView(void *obj, PetscViewer v) {
  PetscViewerASCIIPrintf(v, "name: PetscPF\n");
  PF pf = (PF)obj;
  PFView(pf, v);
  return 0;
}
PetscErrorCode PfR1MatDestroy(void *obj) {
  PF pf = (PF)obj;
  PFDestroy(&pf);
  return 0;
}

PetscErrorCode EER2MatBSS(void *obj, BSS bss, Mat M) {
  PetscErrorCode ierr;
  OpCxtEE* cxt = (OpCxtEE*)obj;
  ierr = BSSEER2Mat(bss, cxt->q, M);
  return 0;
}
PetscErrorCode EER2MatDVR(void *obj, DVR a, Mat M) {
  PetscErrorCode ierr;
  OpCxtEE* cxt = (OpCxtEE*)obj;
  ierr = DVREER2Mat(a, cxt->q, M);
  return 0;
}
PetscErrorCode EER2MatFD(void *obj, FD a, Mat M) {
  PetscErrorCode ierr;
  OpCxtEE* cxt = (OpCxtEE*)obj;
  ierr = FDEER2Mat(a, cxt->q, M);
  return 0;
}
PetscErrorCode EER2MatView(void *obj, PetscViewer v) {
  OpCxtEE* ee = (OpCxtEE*)obj;
  PetscViewerASCIIPrintf(v, "name: e-e interaction\n");
  PetscViewerASCIIPrintf(v, "q: %d\n", ee->q);
  return 0;
}
PetscErrorCode EER2MatDestroy(void *obj) {
  OpCxtEE* ee = (OpCxtEE*)obj;
  PetscFree(ee);
  return 0;
}

PetscErrorCode OpR1Create(MPI_Comm comm, OpR1* p_self) {
  OpR1 self;
  PetscMalloc1(1, &self);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode OpR1Destroy(OpR1* p_self) {

  PetscErrorCode ierr;
  OpR1 self = *p_self;
  if(self->Destroy) {
    ierr = self->Destroy(self->obj); CHKERRQ(ierr);
  }

  ierr = PetscFree(self); CHKERRQ(ierr);

  return 0;  
}

PetscErrorCode OpR1View(OpR1 self, PetscViewer v) {

  PetscViewerASCIIPrintf(v, ">>>> OP(R1) >>>>\n");
  if(self->View)
    self->View(self->obj, v);
  PetscViewerASCIIPrintf(v, "<<<< OP(R1) <<<<\n");  
  return 0;
}

PetscErrorCode OpR1SetD2(OpR1 self) {
  self->obj = NULL;
  self->CalcBSS = D2R1MatBSS;
  self->CalcDVR = D2R1MatDVR;
  self->CalcFD  = D2R1MatFD;
  self->View = D2R1MatView;
  self->Destroy = NULL;
  return 0;
}
PetscErrorCode OpR1SetPf(OpR1 self) {
  self->obj = NULL;
  self->CalcBSS = PfR1MatBSS;
  self->CalcDVR = PfR1MatDVR;
  self->CalcFD  = PfR1MatFD;
  self->View    = PfR1MatView;
  self->Destroy = NULL;
  return 0;
}

PetscErrorCode OpR1Set
