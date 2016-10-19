#include <rescol/bps.h>

// ---- Basic ----
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
PetscErrorCode BPSCopy(BPS self, BPS other) {

  strcpy(other->type, self->type);
  int n = self->num_zs;
  other->num_zs = n;
  PetscMalloc1(n, &other->zs);
  for(int i = 0; i < n; i++) {
    other->zs[i] = self->zs[i];
  }
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
  PetscBool iascii, isbinary, isdraw;

  ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);

  if(v == NULL)
    PetscViewerASCIIGetStdout(self->comm, &v);

  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);  

  if(iascii) {
    PetscViewerASCIIPrintf(v, "BPS object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "type:   %s\n", self->type);
    PetscViewerASCIIPrintf(v, "num_zs: %d\n", self->num_zs);
    PetscViewerASCIIPrintf(v, "zmax:   %f\n", self->zs[self->num_zs-1]);  
    if(self->num_zs > 10) {
      PetscReal* zs = self->zs;
      int n = self->num_zs;
      PetscViewerASCIIPrintf(v, "zs = %f, %f, %f, ..., %f, %f\n",
			     zs[0], zs[1], zs[2], zs[n-1-1], zs[n-1]);
    }
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;

}

// ---- Setter ----
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


// ---- Getter ----
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
PetscErrorCode BPSGetEdge(BPS self, int iele, PetscReal *z0, PetscReal *z1) {

  /*
  PetscErrorCode ierr;
  //ierr = BPSCheckPreallocated(self); CHKERRQ(ierr);
  int num_ele;
  ierr = BPSGetNumEle(self, &num_ele); CHKERRQ(ierr);
  if(iele < 0 || num_ele <= iele ) {
    SETERRQ(self->comm, 1, "iele out of range");
  }
  */

  *z0 = self->zs[iele];
  *z1 = self->zs[iele+1];
  return 0;
}
PetscErrorCode BPSInElementQ(BPS self, int iele, PetscReal x, PetscBool *in_q) {
  // is x in iele element meaning of closed interval.
  double eps = 0.0000001;
  PetscReal z0, z1;
  BPSGetEdge(self, iele, &z0, &z1);
  *in_q = (z0-eps < x && x < z1+eps);
  return 0;
}
