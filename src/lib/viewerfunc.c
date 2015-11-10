#include <rescol/viewerfunc.h>

PetscErrorCode ViewerFuncCreate(ViewerFunc *p_self, MPI_Comm comm) {
  PetscErrorCode ierr;
  ViewerFunc self;
  ierr = PetscNew(&self); CHKERRQ(ierr);

  self->base = NULL;
  
  self->comm = comm;
  self->num = 0;
  self->xs = NULL;
  self->active_base = PETSC_FALSE;
  self->active_range= PETSC_FALSE;
  *p_self = self;

  return 0;
}
PetscErrorCode ViewerFuncDestroy(ViewerFunc *p_self) {
  PetscErrorCode ierr;
  ViewerFunc self = *p_self;

  ierr = PetscViewerDestroy(&self->base); CHKERRQ(ierr);
  ierr = PetscFree(self->xs); CHKERRQ(ierr);

  //ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode ViewerFuncView(ViewerFunc self, PetscViewer viewer) {
 
  PetscViewerType type;
  PetscViewerGetType(viewer, &type);

  if(strcmp(type, "ascii") == 0) {

    FILE *fp;
    PetscViewerASCIIGetPointer(viewer, &fp);
    PetscFPrintf(self->comm, fp, ">>>> ViewerFunc >>>>\n");
    PetscFPrintf(self->comm, fp, "num: %d\n", self->num);
    PetscFPrintf(self->comm, fp, "x[0]: %d\n", self->xs[0]);
    PetscFPrintf(self->comm, fp, "x[num-1]: %d\n", self->xs[self->num-1]);
    PetscFPrintf(self->comm, fp, "<<<< ViewerFunc <<<<\n");

  } else {
    SETERRQ(self->comm, 1, "unsupported type");
  }
  
 
  return 0;
}

PetscErrorCode ViewerFuncSetBase(ViewerFunc self, PetscViewer base) {
  
  self->base = base;
  self->active_base = PETSC_TRUE;

  return 0;
}
PetscErrorCode ViewerFuncSetRange(ViewerFunc self, int num, PetscReal xmax) {

  printf("%d, %f\n", num, xmax);

  self->num = num;
  PetscReal h = xmax/num;  
  PetscMalloc1(num, &self->xs);
  for(int i = 0; i < num; i++) {
    self->xs[i] = (i+1) * h;
  }

  self->active_range = PETSC_TRUE;
  return 0;

}
PetscErrorCode ViewerFuncSetFromOptions(ViewerFunc self) {

  PetscErrorCode ierr;
  PetscBool find, find_xmax, find_num;
  PetscViewer viewer;
  PetscViewerFormat format;

  ierr = PetscOptionsGetViewer(self->comm, NULL, "-viewerfunc_view", 
			       &viewer, &format, &find); CHKERRQ(ierr);
  if(find) {
    ierr = ViewerFuncSetBase(self, viewer); CHKERRQ(ierr);
  }

  PetscInt num;
  ierr = PetscOptionsGetInt(NULL, "-viewerfunc_num",
			    &num, &find_num); CHKERRQ(ierr);

  PetscReal xmax;
  ierr = PetscOptionsGetReal(NULL, "-viewerfunc_xmax",
			     &xmax, &find_xmax); CHKERRQ(ierr);  

  if(find_xmax && find_num) {
    ierr = ViewerFuncSetRange(self ,num, xmax); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode ViewerFuncGetXs(ViewerFunc self, int *num, PetscReal **xs) {
  *num = self->num;
  *xs = self->xs;
  return 0;
}
PetscBool ViewerFuncIsActive(ViewerFunc self) {
  return self->active_range && self->active_base;
}
PetscErrorCode ViewerFuncCheckAcrive(ViewerFunc self) {
  if(ViewerFuncIsActive(self))
    SETERRQ(self->comm, 1, "ViewerFunc is not active");
  return 0;
}
