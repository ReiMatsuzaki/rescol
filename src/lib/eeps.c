#include <rescol/eeps.h>

PetscErrorCode EEPSCreate(MPI_Comm comm, EEPS *p_self) {

  EEPS self;
  PetscNew(&self);

  EPSCreate(comm, &self->eps);
  self->S = NULL;
  self->viewer_values = NULL;
  
  *p_self = self;

  return 0;
}
PetscErrorCode EEPSDestroy(EEPS *p_self) {
  PetscErrorCode ierr;
  EEPS self = *p_self;
  ierr = EPSDestroy(&self->eps); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&self->viewer_values); CHKERRQ(ierr);

  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode EEPSSetOperators(EEPS self, Mat H, Mat S ) {

  PetscErrorCode ierr;
  ierr = EPSSetOperators(self->eps, H, S); CHKERRQ(ierr);
  
  self->S = S;

#if defined(PETSC_USE_COMPLEX)
  if(S == NULL) 
    EPSSetProblemType(self->eps, EPS_NHEP);
  else
    EPSSetProblemType(self->eps, EPS_GNHEP);
#else
  if(S == NULL) 
    EPSSetProblemType(self->eps, EPS_HEP);
  else
    EPSSetProblemType(self->eps, EPS_GHEP);
#endif
  return 0;
}
PetscErrorCode EEPSSetTarget(EEPS self, PetscScalar target) {

  EPSSetWhichEigenpairs(self->eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(self->eps, target);

  return 0;
}
PetscErrorCode EEPSReadInitSpaceOne(EEPS self, PetscViewer viewer) {

  PetscErrorCode ierr;

  Vec x[1]; 
  ierr = VecCreate(PetscObjectComm((PetscObject)self->eps), &x[0]);  CHKERRQ(ierr);

  /*
  printf("EEPSREadInitSpaceOne\n");
  PetscViewer alt_v;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, "evec.dat", FILE_MODE_READ, &alt_v); CHKERRQ(ierr);
  ierr = PetscViewerView(viewer, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = VecLoad(x[0], alt_v); CHKERRQ(ierr);
  */
  ierr = VecLoad(x[0], viewer); CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(self->eps, 1, x); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode EEPSSetInitSpaceFromOther(EEPS self, int n, Mat H, EPS other) {
  
  PetscErrorCode ierr;
  int nconv; EPSGetConverged(other, &nconv);

  if(n == PETSC_DEFAULT)
    n = nconv;

  int n_init = nconv > n ? n : nconv;        // get smaller number

  Vec *xs; PetscMalloc1(n_init, &xs);
  for(int i = 0; i < n_init; i++) {
    ierr = MatCreateVecs(H, &xs[i], NULL);
    ierr = EPSGetEigenpair(other, i, NULL, NULL, xs[i], NULL); CHKERRQ(ierr);
  }  
  ierr = EPSSetInitialSpace(self->eps, n_init, xs); CHKERRQ(ierr);
  for(int i = 0; i < n_init; i++) {
    ierr = VecDestroy(&xs[i]); CHKERRQ(ierr);
    PetscFree(xs);
  }
  return 0;
}
PetscErrorCode EEPSSetFromOptions(EEPS self) {

  PetscErrorCode ierr;
  PetscBool find;
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self->eps, &comm);

  ierr = EPSSetFromOptions(self->eps); CHKERRQ(ierr);

  /* ** I failed to read Vec binary file using PetscOptionsGetViewer written below. **
  PetscViewer init_viewer;
  PetscViewerFormat format;
  ierr = PetscOptionsGetViewer(comm, NULL, "-eeps_view_init", 
			       &init_viewer, &format, &find); CHKERRQ(ierr);
  if(find) {
    ierr = EEPSReadInitSpaceOne(self, init_viewer); CHKERRQ(ierr);
  }
  */ 
  char path[256];
  ierr = PetscOptionsGetString(NULL, "-eeps_path_init", path, 256, &find); CHKERRQ(ierr);
  if(find) {
    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm, path, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = EEPSReadInitSpaceOne(self, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  }
  PetscViewerFormat format;
  ierr = PetscOptionsGetViewer(comm, NULL, "-eeps_view_values", 
			       &self->viewer_values, &format, &find); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode EEPSViewValues(EEPS self) {

  PetscViewerType type;
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self->eps, &comm);
  PetscViewerGetType(self->viewer_values, &type);

  if(strcmp(type, "ascii") != 0) 
    SETERRQ(comm, 1, "only ascii is supported");

  FILE *fp;
  PetscViewerASCIIGetPointer(self->viewer_values, &fp);

  int nconv;
  EPSGetConverged(self->eps, &nconv);
  for(int i = 0 ; i < nconv; i++ ) {
    PetscScalar kr, ki;
    PetscReal re, im;
    PetscReal error;
    EPSGetEigenpair(self->eps, i, &kr, &ki, NULL, NULL);
    EPSComputeError(self->eps, i, EPS_ERROR_RELATIVE, &error);
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif
    PetscFPrintf(comm, fp, "%d %20.16g %20.16g %g\n", i, re, im, error);
		 
  }

  return 0;  
}
PetscErrorCode EEPSSolve(EEPS self) {
  PetscErrorCode ierr;
  ierr = EPSSolve(self->eps); CHKERRQ(ierr);

  if(self->viewer_values != NULL)
    EEPSViewValues(self);
    
  return 0;
}
PetscErrorCode EEPSGetEigenvector(EEPS self, int i, Vec vec) {

  EPSGetEigenpair(self->eps, i, NULL, NULL, vec, NULL);
  
  int n;
  VecGetSize(vec, &n);

#if defined(PETSC_USE_COMPLEX)
  for(int k = 0; k < n; k++) {
    PetscScalar v[1]; PetscInt idx[1] = {k};
    VecGetValues(vec, 1, idx, v);
    if(cabs(v[0]) > 0.00000001) {
      PetscScalar scale = v[0] / cabs(v[0]);
      VecScale(vec, 1.0/scale);
      goto end;
    }
  }
 end:
#endif

  if(self->S == NULL) {
    PetscScalar x;
    VecTDot(vec, vec, &x);
    VecScale(vec, 1.0/sqrt(x));
  } else {
    Vec Sx; MatCreateVecs(self->S, &Sx, NULL);
    PetscScalar scale;
    MatMult(self->S, vec, Sx); VecTDot(vec, Sx, &scale);
    VecScale(vec, 1.0/sqrt(scale));
    VecDestroy(&Sx);
  }

  return 0;
}
