#include <rescol/fd.h>

// ---- Basic Methods ----
PetscErrorCode FDCreate(MPI_Comm comm, FD *p_self) {
  FD self;
  PetscErrorCode ierr;
  ierr = PetscMalloc1(1, &self); CHKERRQ(ierr);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode FDDestroy(FD *p_self) {
  
  PetscErrorCode ierr;
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;

}

PetscErrorCode FDView(FD self, PetscViewer v) {

  PetscViewerType type;
  PetscViewerGetType(v, &type);

  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "unsupported type");
  
  PetscViewerASCIIPrintf(v, ">>>> Finite Difference >>>>\n");
  PetscViewerASCIIPrintf(v, "h  : %f\n", self->h);
  PetscViewerASCIIPrintf(v, "num: %d\n", self->num);
  PetscViewerASCIIPrintf(v, "xmax: %f\n", self->num*self->h);
  PetscViewerASCIIPrintf(v, "<<<< Finite Difference <<<<\n");
  return 0;

}

// ---- Accessor ----
PetscErrorCode FDSetMesh(FD self, int num_xs, PetscReal xmax) {

  self->h = xmax/num_xs;
  self->num = num_xs;

  return 0;
}
PetscErrorCode FDSetFromOptions(FD self) {

  PetscBool find;
  PetscReal xmax = 30.0;  
  PetscInt num = 30;
  PetscErrorCode ierr;

  ierr = PetscOptionsGetReal(NULL, "-fd_xmax", &xmax, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-fd_num", &num, &find); CHKERRQ(ierr);
  
  ierr = FDSetMesh(self, num, xmax); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode FDGetSize(FD self, int *n) {
  *n = self->num;
  return 0;
}

// ---- Matrix/Vector ----
PetscErrorCode FDCreateR1Mat(FD self, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n; FDGetSize(self, &n);
  ierr = MatCreate(self->comm, M); CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
  ierr = MatSetUp(*M); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode FDCreateR2Mat(FD self, Mat *M) {

  PetscErrorCode ierr;
  PetscInt n = self->num;
  ierr = MatCreate(self->comm, M); CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n*n, n*n); CHKERRQ(ierr);
  ierr = MatSetUp(*M); CHKERRQ(ierr);
  return 0;  

}
PetscErrorCode FDCreateR1Vec(FD self, Vec *v) {

  VecCreate(self->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, self->num);
  VecSetUp(*v);
  return 0;

}

PetscErrorCode FDSR1Mat(FD self, Mat M) {

  PetscErrorCode ierr;
  for(int i = 0; i < self->num; i++) {
    ierr = MatSetValue(M, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode FDD2R1Mat(FD self, Mat M) {

  PetscErrorCode ierr;
  PetscInt n = self->num;
  PetscScalar h = self->h;
  PetscScalar hh = h*h;

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      if(i-1 == j || i == j-1)
	ierr = MatSetValue(M, i, j, 1.0/hh, INSERT_VALUES);      
      if(i == j)
	ierr = MatSetValue(M, i, j, -2.0/hh, INSERT_VALUES);
    }
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
  
}
PetscErrorCode FDR2invR1Mat(FD self, Mat M) {
  
  PetscErrorCode ierr;

  for(int i = 0; i < self->num; i++) {
    PetscScalar x = (i+1) * self->h;
    ierr = MatSetValue(M, i, i, 1.0/(x*x), INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode FDENR1Mat(FD self, int q, PetscScalar a, Mat M) {
  
  PetscErrorCode ierr;

  for(int i = 0; i < self->num; i++) {
    PetscReal x = (i+1) * self->h;
    PetscReal y;
    PartialCoulomb(q, a, x, &y);
    ierr = MatSetValue(M, i, i, y, INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
  
}
PetscErrorCode FDEER2Mat(FD self, int q, Mat M) {

  PetscErrorCode ierr;
  PetscInt n = self->num;
  PetscScalar h = self->h;

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal xi = (i+1) * h;
      PetscReal xj = (j+1) * h;
      PetscReal y;
      PartialCoulomb(q, xi, xj, &y);
      PetscInt idx = i*n+j;
      ierr = MatSetValue(M, idx, idx, y, INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode FDGuessHEig(FD self, int n, int l, PetscScalar z, Vec v) {
  
  if(n != 1 || l != 0) 
    SETERRQ(self->comm, 1, "only (n=1,l=0) is supported now.");

  for(int i = 0; i < self->num; i++) {
    PetscScalar x = (i+1)*self->h;
    PetscScalar y = 2.0*pow(z, 1.5)*x*exp(-z*x);
    VecSetValue(v, i, y, INSERT_VALUES);
  }
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  return 0;

}
