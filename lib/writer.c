#include <rescol/writer.h>

PetscErrorCode WFWriterCreate(WFWriter *p_self, MPI_Comm comm) {
  WFWriter self;
  PetscNew(&self);
  self->comm = comm;
  self->active = PETSC_FALSE;
  self->num = 0;
  self->xmax = 0.0;
  self->fp = NULL;
  *p_self = self;
  return 0;
}
PetscErrorCode WFWriterSet(WFWriter self, PetscInt num, PetscReal xmax) {
  self->active = PETSC_TRUE;
  self->num = num; self->xmax=xmax; 
  return 0;
}
PetscErrorCode WFWriterSetPath(WFWriter self, char path[]) {
  PetscErrorCode ierr;
  ierr = PetscFOpen(self->comm, path, "w", &self->fp); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode WFWriterSetFromOptions(WFWriter self) {

  PetscBool find_num, find_xmax, find_path;
  PetscInt num;
  PetscReal xmax;
  char path[256];

  PetscOptionsGetInt(NULL, "-wfwriter_num", &num, &find_num);
  PetscOptionsGetReal(NULL, "-wfwriter_xmax", &xmax, &find_xmax);
  PetscOptionsGetString(NULL, "-wfwriter_view", path, 256, &find_path);

  if(!find_num || !find_xmax || !find_path) {
    return 0;
  }

  PetscErrorCode ierr;
  ierr = WFWriterSet(self, num, xmax); CHKERRQ(ierr);
  ierr = WFWriterSetPath(self, path); CHKERRQ(ierr);

  return 0;
}
/*
PetscErrorCode WFWriterCreateFromOptions(WFWriter *p_self, MPI_Comm comm) {
  PetscErrorCode ierr;
  ierr = WFWriterCreate(p_self, comm); CHKERRQ(ierr);
  ierr = WFWriterSetFromOptions(*p_self); CHKERRQ(ierr);
  return 0;
}
*/
PetscErrorCode WFWriterDestroy(WFWriter *p_self) {
  PetscErrorCode ierr;
  WFWriter self = *p_self;
  if(self->fp != NULL) {
    ierr = PetscFClose(self->comm, self->fp); CHKERRQ(ierr);
  }
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
  PetscPrintf(self->comm, "f: %f\n", self->xmax);  
  PetscPrintf(self->comm, "<<<< WFWriter <<<<\n");

  return 0;
}
PetscBool WFWriterIsActive(WFWriter self) {
  return self->active;
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

