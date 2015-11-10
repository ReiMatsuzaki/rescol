#include <rescol/eeps.h>

PetscErrorCode EEPSCreate(EEPS *p_self, MPI_Comm comm) {

  EEPS self;
  PetscNew(&self);

  self->comm = comm;
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

  EPSSetOperators(self->eps, H, S);
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

  Vec x[1]; VecCreate(self->comm, &x[0]); 
  ierr = VecLoad(x[0], viewer); CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(self->eps, 1, x); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode EEPSSetFromOptions(EEPS self) {

  PetscErrorCode ierr;
  PetscBool find;

  ierr = EPSSetFromOptions(self->eps); CHKERRQ(ierr);
 
  PetscViewer init_viewer;
  PetscViewerFormat format;
  ierr = PetscOptionsGetViewer(self->comm, NULL, "-eeps_view_init", 
			       &init_viewer, &format, &find); CHKERRQ(ierr);
  if(find) {
    ierr = EEPSReadInitSpaceOne(self, init_viewer); CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetViewer(self->comm, NULL, "-eeps_view_values", 
			       &self->viewer_values, &format, &find); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode EEPSViewValues(EEPS self) {

  PetscViewerType type;
  PetscViewerGetType(self->viewer_values, &type);

  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "only ascii is supported");

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
    im = PetscRealPart(kr);
#else
    re = kr;
    im = ki;
#endif
    PetscFPrintf(self->comm, fp, "%d %20.16g %20.16g %g\n", i, re, im, error);
		 
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
  }

  return 0;
}