#include <rescol/writer.h>

PetscErrorCode WFWriterCreate(WFWriter *self, MPI_Comm comm) {
  WFWriter _self;
  PetscNew(&_self);
  _self->comm = comm;
  *self = _self;
  return 0;
}
PetscErrorCode WFWriterSet(WFWriter self, PetscInt num, PetscReal xmax) {
  self->num = num; self->xmax=xmax;
  return 0;
}
PetscErrorCode WFWriterCreateFromOptions(WFWriter *p_self, MPI_Comm comm) {
  PetscBool find_num, find_xmax;
  PetscInt num;
  PetscReal xmax;

  PetscOptionsGetInt(NULL, "-wfwriter_num", &num, &find_num);
  PetscOptionsGetReal(NULL, "-wfwriter_xmax", &xmax, &find_xmax);

  if(!find_num || !find_xmax) {
    *p_self = NULL;
    return 0;
  }

  WFWriterCreate(p_self, comm);
  WFWriterSet(*p_self, num, xmax);
  return 0;
}
PetscErrorCode WFWriterDestroy(WFWriter *p_self) {
  PetscFree(*p_self);
  return 0;
}

PetscErrorCode WFWriterCheckState(WFWriter self) {
  if(self==NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "self is null");
  return 0;
}
PetscErrorCode WFWriterView(WFWriter self) {

  PetscErrorCode ierr;
  ierr = WFWriterCheckState(self); CHKERRQ(ierr);

  PetscPrintf(self->comm, ">>>> WFWriter >>>>\n");
  PetscPrintf(self->comm, "num: %d\n", self->num);
  PetscPrintf(self->comm, "xmax: %f\n", self->xmax);  
  PetscPrintf(self->comm, "<<<< WFWriter <<<<\n");

  return 0;
}

PetscErrorCode WFWriterWrite(WFWriter self, FILE *fp, FEMInf fem, Vec c) {

  PetscErrorCode ierr;
  ierr = WFWriterCheckState(self); CHKERRQ(ierr);
  MPI_Comm comm = self->comm;

  if(fp == NULL)
    SETERRQ(comm, 1, "fp is null");

  PetscReal h = self->xmax/self->num;

  for(int i = 0; i <= self->num; i++) {
    PetscReal x = i * h;
    PetscScalar y;
    ierr = FEMInfPsi(fem, x, c, &y); CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    PetscReal yr = PetscRealPart(y);
    PetscReal yi = PetscImaginaryPart(y);
    PetscFPrintf(comm, fp, "%f %f %f\n", x, yr, yi);
#else
    PetscFPrintf(comm, fp, "%f %f\n", x, y);
#endif
    
  }
  
  return 0;
}
PetscErrorCode WFWriterWriteFile(WFWriter self, char *fn, FEMInf fem, Vec c) {
  PetscErrorCode ierr;
  FILE *fp;
  ierr = PetscFOpen(self->comm, fn, "w", &fp); CHKERRQ(ierr);
  ierr = WFWriterWrite(self, fp, fem, c); CHKERRQ(ierr);
  ierr = PetscFClose(self->comm, fp);   CHKERRQ(ierr);
  return 0;
}

PetscErrorCode WFWriterWriteOCE1(WFWriter self, FILE *fp, OCE1 oce, Vec c) {

  PetscErrorCode ierr;
  ierr = WFWriterCheckState(self); CHKERRQ(ierr);

  int n_y, n_r; OCE1GetSizes(oce, &n_r, &n_y);
  if(n_y != 1)
    SETERRQ(self->comm, 1, "now only n_y == 1 is supported");

  WFWriterWrite(self, fp, oce->fem, c);
  return 0;
}
PetscErrorCode WFWriterWriteFileOCE1(WFWriter self, char *fn, OCE1 oce, Vec c) {
  PetscErrorCode ierr;
  FILE *fp;
  ierr = PetscFOpen(self->comm, fn, "w", &fp); CHKERRQ(ierr);
  ierr = WFWriterWriteOCE1(self, fp, oce, c); CHKERRQ(ierr);
  ierr = PetscFClose(self->comm, fp);   CHKERRQ(ierr);
  return 0;

}

PetscErrorCode WFWriterWritePOT(WFWriter self, FILE *fp, POT pot, int J, PetscReal mu) {

  PetscErrorCode ierr;
  MPI_Comm comm = self->comm;
  ierr = WFWriterCheckState(self); CHKERRQ(ierr);

  if(fp == NULL)
    SETERRQ(comm, 1, "fp is null");

  PetscReal h = self->xmax/self->num;

  for(int i = 0; i <= self->num; i++) {
    PetscReal x = i * h;
    PetscScalar y; POTCalc(pot, x, &y);
    if(i != 0)
      y += J*(J+1)/(2.0*mu*x*x);
#if defined(PETSC_USE_COMPLEX)
    PetscReal yr = PetscRealPart(y);
    PetscReal yi = PetscImaginaryPart(y);
    PetscFPrintf(comm, fp, "%f %f %f\n", x, yr, yi);
#else
    PetscFPrintf(comm, fp, "%f %f\n", x, y);
#endif
    
  }  

  return 0;
}
PetscErrorCode WFWriterWriteFilePOT(WFWriter self, char *fn, POT pot, int J, PetscReal mu) {

  PetscErrorCode ierr;
  FILE *fp;
  ierr = PetscFOpen(self->comm, fn, "w", &fp); CHKERRQ(ierr);
  ierr = WFWriterWritePOT(self, fp, pot, J, mu); CHKERRQ(ierr);
  ierr = PetscFClose(self->comm, fp);   CHKERRQ(ierr);
  return 0;

}

