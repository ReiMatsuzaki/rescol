#include <rescol/bps.h>

// ------- Basic ----------
PetscErrorCode BPSCreate(BPS *bps, MPI_Comm comm) {
  //  PetscErrorCode ierr;
  BPS _bps;
  PetscNew(&_bps);
  *bps = NULL;

  strcpy(_bps->type, "unknows");
  _bps->num_zs = 0;
  _bps->zs = NULL;
  _bps->comm = comm;

  *bps = _bps;
  return 0;
}

PetscErrorCode BPSSetExp(BPS this, PetscScalar zmax, PetscInt num_zs, PetscScalar gamma) {

  if(this == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "BPS object is null");

  PetscErrorCode ierr;
  strcpy(this->type, "exp");
  this->num_zs = num_zs;
  ierr = PetscMalloc1(num_zs, &this->zs); CHKERRQ(ierr);
  for(int i = 0; i < num_zs; i++)
    this->zs[i] = zmax * (exp(gamma*i/(num_zs-1)) - 1.0) / (exp(gamma) - 1.0);

  return 0;  

}

PetscErrorCode BPSSetLine(BPS this, PetscScalar zmax, PetscInt num_zs) {

  if(this == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "BPS object is null");

  PetscErrorCode ierr;
  strcpy(this->type, "line");
  this->num_zs = num_zs;
  ierr = PetscMalloc1(num_zs, &this->zs);
  for(int i = 0; i < num_zs; i++)
    this->zs[i] = i * zmax / (num_zs - 1);

  return 0;

}
PetscErrorCode BPSSetFromOptions(BPS this) {

  PetscBool find;
  PetscReal zmax = 20.0;
  PetscInt num_zs = 3;
  char type[10] = "line";
  PetscErrorCode ierr;

  ierr = PetscOptionsGetInt(NULL, "-bps_num_zs", &num_zs, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, "-bps_zmax", &zmax, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-bps_type", type, 10, &find); CHKERRQ(ierr);

  if(strcmp(type, "line") == 0 ) {
    ierr = BPSSetLine(this, zmax, num_zs); CHKERRQ(ierr); 
  } else if(strcmp(type, "exp") == 0) {
    ierr = BPSSetExp(this, zmax, num_zs, 5.0); CHKERRQ(ierr);
  } else
    SETERRQ(this->comm, 1, "-bps_type must be line or exp");

  return 0;

}
PetscErrorCode BPSDestroy(BPS *bps) {
  PetscErrorCode ierr;
  ierr = PetscFree((*bps)->zs); CHKERRQ(ierr);
  ierr = PetscFree(*bps); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BPSCheckPreallocated(BPS this) {

  if(this == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "bps object is null");

  if(this->zs == NULL)
    SETERRQ(this->comm, 1, "BPS object is not allocated. call BPSSet*** function before using this object");

  return 0;

}
PetscErrorCode BPSFPrintf(BPS this, FILE *file, int lvl) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(this); CHKERRQ(ierr);
  
  if(lvl != 0) 
    SETERRQ(this->comm, 1, "now only lvl=0 is supported");

  PetscFPrintf(this->comm, file, "==== Begin BPS ====\n");
  PetscFPrintf(this->comm, file, "type: %s\n", this->type);
  PetscFPrintf(this->comm, file, "# of points: %d\n", this->num_zs);
  PetscFPrintf(this->comm, file, "zmax: %f\n", this->zs[this->num_zs-1]);  
  PetscFPrintf(this->comm, file, "==== End BPS ====\n");
  return 0;

}

// ------ Getter ---------
PetscErrorCode BPSGetZs(BPS this, PetscScalar **zs, PetscInt *num_zs) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(this); CHKERRQ(ierr);
  
  if(num_zs != NULL)
    *num_zs = this->num_zs;
  
  ierr = PetscMalloc1(this->num_zs, zs);
  for(int i = 0; i < this->num_zs; i++)
    (*zs)[i] = this->zs[i];

  return 0;
}
PetscErrorCode BPSGetNumEle(BPS this, PetscInt *num_ele) {

  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(this); CHKERRQ(ierr);

  *num_ele = this->num_zs - 1;

  return 0;
}
PetscErrorCode BPSGetZMax(BPS this, PetscScalar *zmax) {
  PetscErrorCode ierr;
  ierr = BPSCheckPreallocated(this); CHKERRQ(ierr);
  *zmax = this->zs[this->num_zs-1];
  return 0;
}
