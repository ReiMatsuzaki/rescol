#include <rescol/bps.h>

// ------- Basic ----------
PetscErrorCode BPSCreate(MPI_Comm comm, BPS *p_self) {
  //  PetscErrorCode ierr;
  BPS self;
  PetscNew(&self);

  strcpy(self->type, "unknows");
  self->num_zs = 0;
  self->zs = NULL;
  self->comm = comm;

  *p_self = self;
  return 0;
}
PetscErrorCode BPSDestroy(BPS *p_self) {
  PetscErrorCode ierr;
  ierr = PetscFree((*p_self)->zs); CHKERRQ(ierr);
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BPSCheckPreallocated(BPS self) {

  if(self == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "bps object is null");

  if(self->zs == NULL)
    SETERRQ(self->comm, 1, "BPS object is not allocated. call BPSSet*** function before using self object");

  return 0;

}
PetscErrorCode BPSView(BPS self, PetscViewer v) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);

  PetscViewerASCIIPrintf(v, ">>>> BPS >>>>\n");
  PetscViewerASCIIPrintf(v, "type: %s\n", self->type);
  PetscViewerASCIIPrintf(v, "num_of_points: %d\n", self->num_zs);
  PetscViewerASCIIPrintf(v, "zmax: %f\n", self->zs[self->num_zs-1]);  
  PetscViewerASCIIPrintf(v, "<<<< BPS <<<<\n");
  return 0;

}


PetscErrorCode BPSSetExp(BPS self, PetscReal zmax, PetscInt num_zs, PetscReal gamma) {
  
  PetscErrorCode ierr;
  strcpy(self->type, "exp");
  self->num_zs = num_zs;
  ierr = PetscMalloc1(num_zs, &self->zs); CHKERRQ(ierr);
  for(int i = 0; i < num_zs; i++)
    self->zs[i] = zmax * (exp(gamma*i/(num_zs-1)) - 1.0) / (exp(gamma) - 1.0);

  return 0;  

}

PetscErrorCode BPSSetLine(BPS self, PetscReal zmax, PetscInt num_zs) {

  if(self == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "BPS object is null");

  PetscErrorCode ierr;
  strcpy(self->type, "line");
  self->num_zs = num_zs;
  ierr = PetscMalloc1(num_zs, &self->zs);
  for(int i = 0; i < num_zs; i++)
    self->zs[i] = i * zmax / (num_zs - 1);

  return 0;

}
PetscErrorCode BPSSetFromOptions(BPS self) {

  PetscBool find;
  PetscReal zmax = 20.0;
  PetscInt num_zs = 3;
  char type[10] = "line";
  PetscErrorCode ierr;

  ierr = PetscOptionsGetInt(NULL, "-bps_num_zs", &num_zs, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, "-bps_zmax", &zmax, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-bps_type", type, 10, &find); CHKERRQ(ierr);

  if(strcmp(type, "line") == 0 ) {
    ierr = BPSSetLine(self, zmax, num_zs); CHKERRQ(ierr); 
  } else if(strcmp(type, "exp") == 0) {
    ierr = BPSSetExp(self, zmax, num_zs, 5.0); CHKERRQ(ierr);
  } else
    SETERRQ(self->comm, 1, "-bps_type must be line or exp");

  return 0;

}


// ------ Getter ---------
PetscErrorCode BPSGetZs(BPS self, PetscReal **zs, PetscInt *num_zs) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);
  
  if(num_zs != NULL)
    *num_zs = self->num_zs;
  
  ierr = PetscMalloc1(self->num_zs, zs);
  for(int i = 0; i < self->num_zs; i++)
    (*zs)[i] = self->zs[i];

  return 0;
}
PetscErrorCode BPSGetNumEle(BPS self, PetscInt *num_ele) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);

  *num_ele = self->num_zs - 1;

  return 0;
}
PetscErrorCode BPSGetZMax(BPS self, PetscReal *zmax) {
  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);
  *zmax = self->zs[self->num_zs-1];
  return 0;
}
