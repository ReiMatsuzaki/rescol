#include "../include/viewerfunc.h"

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
  strcpy(self->opt_prefix, "");
  strcpy(self->opt_base, "viewerfunc");
  *p_self = self;

  return 0;
}
PetscErrorCode ViewerFuncDestroy(ViewerFunc *p_self) {
  PetscErrorCode ierr=0;
  ViewerFunc self = *p_self;

  if(self->base != NULL)
    ierr = PetscViewerDestroy(&self->base); CHKERRQ(ierr);

  if(self->xs != NULL)
    ierr = PetscFree(self->xs); CHKERRQ(ierr);

  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscViewer ViewerFuncGetBase(ViewerFunc self) {
  return self->base;
}
PetscErrorCode ViewerFuncView(ViewerFunc self, PetscViewer v) {

  if(!ViewerFuncIsActive(self)) {
    SETERRQ(self->comm, 1, "ViewerFunc object is not setup");
  }

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

    PetscBool is_active = ViewerFuncIsActive(self);
    PetscViewerASCIIPrintf(v, "is_active: %s\n",
			   is_active ? "Yes" : "No");

    if(is_active) {
      PetscViewerASCIIPrintf(v, "num: %d\n", self->num);
      PetscReal *xs = self->xs;
      int n = self->num;
      if(self->num > 5) {
	PetscViewerASCIIPrintf(v, "xs = %f, %f, %f, ..., %f, %f\n",
			       xs[0], xs[1], xs[2], xs[n-2], xs[n-1]);
      } else {
	PetscViewerASCIIPrintf(v, "xs = ");
	for(int i = 0; i < n; i++) {
	  PetscViewerASCIIPrintf(v, "%f", xs[i]);
	  if(i != n-1)
	    PetscViewerASCIIPrintf(v, ", ");
	}
      }
    }
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

  PetscErrorCode ierr;
  self->num = num;
  PetscReal h = xmax/(num-1);

  if(self->xs != NULL) {
    ierr = PetscFree(self->xs); CHKERRQ(ierr);
  }
  
  PetscMalloc1(num, &self->xs);
  for(int i = 0; i < num; i++) {
    self->xs[i] = i * h;
  }

  self->active_range = PETSC_TRUE;
  return 0;

}
PetscErrorCode ViewerFuncSetOptionsPrefix(ViewerFunc self, const char prefix[]) {
  strcpy(self->opt_prefix, prefix);
  return 0;
}
PetscErrorCode ViewerFuncSetFromOptions(ViewerFunc self) {

  PetscErrorCode ierr;

  char opt_path[100];
  char opt_num[100];
  char opt_xmax[100];

  sprintf(opt_path, "-%s%s_path",
	  self->opt_prefix, self->opt_base);
  sprintf(opt_num,  "-%s%s_num",
	  self->opt_prefix, self->opt_base);
  sprintf(opt_xmax, "-%s%s_xmax",
	  self->opt_prefix, self->opt_base);

  PetscBool find_path, find_xmax, find_num;
  char path[100];
  PetscInt  num;
  PetscReal xmax;
  ierr = PetscOptionsGetString(NULL, NULL,opt_path, path,  100, &find_path); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL,   opt_num,  &num, &find_num); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL,  opt_xmax, &xmax, &find_xmax); CHKERRQ(ierr);

  if(find_path) {
    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(self->comm, path, &viewer); CHKERRQ(ierr);
    ierr = ViewerFuncSetBase(self, viewer); CHKERRQ(ierr);
  }

  if(find_num && find_xmax) {
    ierr = ViewerFuncSetRange(self, num, xmax); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode ViewerFuncGetRange(ViewerFunc self, int *num, PetscReal **xs) {

  if(!ViewerFuncIsActive(self))
    SETERRQ(self->comm, 1, "ViewerFunc object is not setup");
  
  *num = self->num;
  *xs = self->xs;
  return 0;
}
PetscBool ViewerFuncIsActive(ViewerFunc self) {
  return self->active_range && self->active_base;
}

