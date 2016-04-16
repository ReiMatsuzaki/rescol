#include <rescol/viewerfunc.h>

PetscErrorCode ViewerFuncCreate(MPI_Comm comm, ViewerFunc *p_self) {
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

  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscViewer ViewerFuncGetBase(ViewerFunc self) {
  return self->base;
}
PetscErrorCode ViewerFuncView(ViewerFunc self, PetscViewer v) {

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

  if(iascii) {
    PetscViewerASCIIPrintf(v, "ViewerFunc object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "num: %d\n", self->num);
    PetscViewerASCIIPrintf(v, "x[0]: %d\n", self->xs[0]);
    PetscViewerASCIIPrintf(v, "x[num-1]: %d\n", self->xs[self->num-1]);
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}

PetscErrorCode ViewerFuncSetBase(ViewerFunc self, PetscViewer base) {
  
  self->base = base;
  self->active_base = PETSC_TRUE;

  return 0;
}
PetscErrorCode ViewerFuncSetRange(ViewerFunc self, int num, PetscReal xmax) {

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
  } else {
    SETERRQ(self->comm, 1, "-viewerfunc_num and -viewerfunc_xmax is necessary");
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
