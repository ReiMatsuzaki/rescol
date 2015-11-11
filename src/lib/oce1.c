#include <rescol/oce1.h>

PetscErrorCode OCE1Create(MPI_Comm comm, OCE1 *p_self) {

  PetscErrorCode ierr;
  OCE1 self;

  ierr = PetscNew(&self); CHKERRQ(ierr);

  self->comm = comm;
  self->mu = 1.0;
  self->fem = NULL;
  self->y1s = NULL;
  self->s_r = NULL;
  self->s_y = NULL;
  
  *p_self = self;
  return 0;
  
}
PetscErrorCode OCE1Destroy(OCE1 *p_self) {

  if((*p_self)->fem)
    FEMInfDestroy(&(*p_self)->fem);

  if((*p_self)->y1s)
    Y1sDestroy(&(*p_self)->y1s);  

  if((*p_self)->s_r)
    MatDestroy(&(*p_self)->s_r);  

  if((*p_self)->s_y)
    MatDestroy(&(*p_self)->s_y);  

  PetscFree(*p_self);
  return 0;

}

PetscErrorCode OCE1View(OCE1 self, PetscViewer v) {

  PetscViewerType type;
  PetscViewerGetType(v, &type);

  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "unsupported type");

  PetscViewerASCIIPrintf(v, ">>>> One Center Expansion 1 >>>>\n");
  PetscViewerASCIIPrintf(v, "mu: %f\n", self->mu);
  FEMInfView(self->fem, v);
  Y1sView(self->y1s, v);
  PetscViewerASCIIPrintf(v, "<<<< One Center Expansion 1 <<<<\n");
  return 0;
}

PetscErrorCode OCE1Set(OCE1 self, FEMInf fem, Y1s y1s) {

  self->fem = fem;
  self->y1s = y1s;
  return 0;

}
PetscErrorCode OCE1SetMu(OCE1 self, PetscReal mu) {
  self->mu = mu;
  return 0;
}
PetscErrorCode OCE1SetFromOptions(OCE1 self) {

  PetscErrorCode ierr;
  PetscBool find;
  PetscReal mu;

  PetscOptionsGetReal(NULL, "-oce1_mu", &mu, &find); 
  if(find)
    OCE1SetMu(self, mu);

  FEMInf fem; FEMInfCreate(self->comm, &fem);
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);

  Y1s y1s; Y1sCreate(self->comm, &y1s);
  ierr = Y1sSetFromOptions(y1s); CHKERRQ(ierr);

  ierr = OCE1Set(self, fem, y1s); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode OCE1GetSizes(OCE1 self, int *n_r, int *n_y) {
  FEMInfGetSize(self->fem, n_r);
  Y1sGetSize(self->y1s, n_y);
  return 0;
}

PetscErrorCode OCE1SetSr(OCE1 self) {
  PetscErrorCode ierr;
  ierr = FEMInfCreateMat(self->fem, 1, &self->s_r);CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(self->fem, self->s_r);CHKERRQ(ierr);
  return 0;
}
PetscErrorCode OCE1SetSy(OCE1 self) {
  PetscErrorCode ierr;
  ierr = Y1sCreateY1Mat(self->y1s, &self->s_y); CHKERRQ(ierr);
  ierr = Y1sSY1Mat(self->y1s, self->s_y);CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode OCE1SetSMat(OCE1 self, Mat *M) {

  PetscErrorCode ierr;

  if(self->s_r == NULL) {
    ierr = OCE1SetSr(self); CHKERRQ(ierr);
  }
  if(self->s_y == NULL) {
    ierr = OCE1SetSy(self); CHKERRQ(ierr);
  }

  ierr = MatSetSynthesizeFast(self->s_r, self->s_y, self->comm, M);
  CHKERRQ(ierr);

  return 0;
}
PetscErrorCode OCE1SetSMatNullable(OCE1 self, Mat *M) {

  PetscBool s_is_id;
  FEMInfGetOverlapIsId(self->fem, &s_is_id);

  if(s_is_id)
    *M = NULL;
  else {
    PetscErrorCode ierr;
    ierr = OCE1SetSMat(self, M); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode OCE1SetTMat(OCE1 self, Mat *M) {

  if(self->s_r == NULL)
    OCE1SetSr(self);
  if(self->s_y == NULL)
    OCE1SetSy(self);

  PetscErrorCode ierr;

  Mat d2_r;
  ierr = FEMInfSetD2R1Mat(self->fem, &d2_r); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(d2_r, self->s_y, self->comm, M); CHKERRQ(ierr);
  ierr = MatScale(*M, -0.5/self->mu); CHKERRQ(ierr);
  ierr = MatDestroy(&d2_r); CHKERRQ(ierr);

  Mat l_r, l_y, L;
  ierr = FEMInfSetR2invR1Mat(self->fem, &l_r); CHKERRQ(ierr);
  ierr = Y1sSetLambdaY1Mat(self->y1s, &l_y); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(l_r, l_y, self->comm, &L);
  ierr = MatAXPY(*M, 0.5/self->mu, L, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&l_r); MatDestroy(&l_y); MatDestroy(&L);

  return 0;
}
PetscErrorCode OCE1PlusPOTMat(OCE1 self, RotSym sym, POT pot, Mat M) {

  if(sym != ROT_SCALAR)
    SETERRQ(self->comm, 1, "now only sym=ROT_SCALAR is supported");

  PetscErrorCode ierr;
  if(self->s_y == NULL)
    OCE1SetSy(self);

  Mat r1, V;
  ierr = FEMInfSetPOTR1Mat(self->fem, pot, &r1); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(r1, self->s_y, self->comm, &V); CHKERRQ(ierr);
  MatDestroy(&r1);

  MatAXPY(M, 1.0, V, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&V);

  return 0;
}
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat M) {

  PetscErrorCode ierr;

  int lmax; Y1sGetMaxL(self->y1s, &lmax);
  int qmax = 2*lmax;
  PetscReal zz = -2.0*z;

  for(int q = 0; q < qmax; q++) {
    Mat pq_y = NULL;
    ierr = Y1sSetPqY1Mat(self->y1s, q, &pq_y); CHKERRQ(ierr);

    if(pq_y) {
      Mat pq_r;
      ierr = FEMInfSetENR1Mat(self->fem, q, a, &pq_r); CHKERRQ(ierr);

      Mat V;
      ierr = MatSetSynthesizeFast(pq_r, pq_y, self->comm, &V); CHKERRQ(ierr);

      ierr = MatAXPY(M, zz, V, DIFFERENT_NONZERO_PATTERN);
      
      MatDestroy(&pq_r); MatDestroy(&pq_y); MatDestroy(&V);
    }
  }
  return 0;
  
}
