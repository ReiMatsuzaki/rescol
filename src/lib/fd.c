#include <rescol/fd.h>

// ---- Basic Methods ----
PetscErrorCode FDCreate(FD *fd, int num_xs, double xmax, MPI_Comm comm) {

  PetscErrorCode ierr;
  FD _fd;
  //  _fd = (FD)malloc(sizeof(struct _p_FD));
  ierr = PetscMalloc1(1, &_fd);
  *fd = NULL;

  _fd->comm = comm;
  _fd->h = xmax/num_xs;
  _fd->num = num_xs;

  *fd = _fd;
  return 0;
}
PetscErrorCode FDCreateFromOptions(FD *fd, MPI_Comm comm) {

  PetscBool find;
  PetscReal xmax = 30.0;  
  PetscInt num = 30;
  PetscErrorCode ierr;

  ierr = PetscOptionsGetReal(NULL, "-fd_xmax", &xmax, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-fd_num", &num, &find); CHKERRQ(ierr);
  
  ierr = FDCreate(fd, num, xmax, comm); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode FDDestroy(FD *fd) {
  
  PetscErrorCode ierr;
  ierr = PetscFree(*fd); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode FDFPrintf(FD this, FILE *file, int lvl) {
  
  if(lvl != 0) 
    SETERRQ(this->comm, 1, "now only lvl=0 is supported");

  PetscFPrintf(this->comm, file, "==== Begin Finite Difference ====\n");
  PetscFPrintf(this->comm, file, "h  : %f\n", this->h);
  PetscFPrintf(this->comm, file, "num: %d\n", this->num);
  PetscFPrintf(this->comm, file, "xmax: %f\n", this->num*this->h);
  PetscFPrintf(this->comm, file, "==== End Finite Difference ====\n");
  return 0;

}

// ---- Accessor ----
PetscErrorCode FDGetSize(FD self, int *n) {
  *n = self->num;
  return 0;
}

// ---- Matrix (private) ----
PetscErrorCode FDInitR1Mat(FD this, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n; FDGetSize(this, &n);
  ierr = MatCreate(this->comm, M); CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*M); CHKERRQ(ierr);
  ierr = MatSetUp(*M); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode FDInitR2Mat(FD this, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n = this->num;
  ierr = MatCreate(this->comm, M); CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n*n, n*n); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*M); CHKERRQ(ierr);
  ierr = MatSetUp(*M); CHKERRQ(ierr);
  return 0;  

}

// ---- Matrix (public) ----
PetscErrorCode FDSetSR1Mat(FD this, Mat *M) {

  PetscErrorCode ierr;
  ierr = FDInitR1Mat(this, M); CHKERRQ(ierr);
  for(int i = 0; i < this->num; i++) {
    ierr = MatSetValue(*M, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode FDSetD2R1Mat(FD this, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n = this->num;
  PetscScalar h = this->h;
  PetscScalar hh = h*h;

  ierr = FDInitR1Mat(this, M); CHKERRQ(ierr);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      if(i-1 == j || i == j-1)
	ierr = MatSetValue(*M, i, j, 1.0/hh, INSERT_VALUES);      
      if(i == j)
	ierr = MatSetValue(*M, i, j, -2.0/hh, INSERT_VALUES);
    }
  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
  
}
PetscErrorCode FDSetR2invR1Mat(FD this, Mat *M) {
  
  PetscErrorCode ierr;

  FDInitR1Mat(this, M);
  for(int i = 0; i < this->num; i++) {
    PetscScalar x = (i+1) * this->h;
    ierr = MatSetValue(*M, i, i, 1.0/(x*x), INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode FDSetENR1Mat(FD this, int q, PetscScalar a, Mat *M) {
  
  PetscErrorCode ierr;

  FDInitR1Mat(this, M);
  for(int i = 0; i < this->num; i++) {
    PetscReal x = (i+1) * this->h;
    PetscReal y;
    PartialCoulomb(q, a, x, &y);
    ierr = MatSetValue(*M, i, i, y, INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
  
}
PetscErrorCode FDSetEER2Mat(FD this, int q, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n = this->num;
  PetscScalar h = this->h;

  FDInitR2Mat(this, M);
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal xi = (i+1) * h;
      PetscReal xj = (j+1) * h;
      PetscReal y;
      PartialCoulomb(q, xi, xj, &y);
      PetscInt idx = i*n+j;
      ierr = MatSetValue(*M, idx, idx, y, INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}

// ---- Vector ----
PetscErrorCode FDGuessHEig(FD this, int n, int l, PetscScalar z, Vec *v) {
  
  if(n != 1 || l != 0) 
    SETERRQ(this->comm, 1, "only (n=1,l=0) is supported now.");

  VecCreate(this->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, this->num);
  VecSetFromOptions(*v);
  VecSetUp(*v);

  for(int i = 0; i < this->num; i++) {
    PetscScalar x = (i+1)*this->h;
    PetscScalar y = 2.0*pow(z, 1.5)*x*exp(-z*x);
    VecSetValue(*v, i, y, INSERT_VALUES);
  }
  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);
  return 0;

}
