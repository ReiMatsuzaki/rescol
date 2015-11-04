#include "oce1.h"

PetscErrorCode OCE1Create(OCE1 *oce1, MPI_Comm comm) {

  PetscErrorCode ierr;
  
  OCE1 _oce1;

  ierr = PetscNew(&_oce1); CHKERRQ(ierr);
  *oce1 = NULL;

  _oce1->comm = comm;
  _oce1->fem = NULL;
  _oce1->y1s = NULL;
  _oce1->s_r = NULL;
  _oce1->s_y = NULL;
  
  *oce1 = _oce1;
  return 0;
  
}
PetscErrorCode OCE1Destroy(OCE1 *oce1) {

  if((*oce1)->fem)
    FEMInfDestroy(&(*oce1)->fem);

  if((*oce1)->y1s)
    Y1sDestroy(&(*oce1)->y1s);  

  if((*oce1)->s_r)
    MatDestroy(&(*oce1)->s_r);  

  if((*oce1)->s_y)
    MatDestroy(&(*oce1)->s_y);  

  PetscFree(*oce1);
  return 0;

}
PetscErrorCode OCE1Set(OCE1 self, FEMInf fem, Y1s y1s) {

  self->fem = fem;
  self->y1s = y1s;
  return 0;

}
PetscErrorCode OCE1CreateFromOptions(OCE1 *oce1, MPI_Comm comm) {

  PetscErrorCode ierr;
  ierr = OCE1Create(oce1, comm); CHKERRQ(ierr);

  FEMInf fem;
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);

  Y1s y1s;
  ierr = Y1sCreateFromOptions(&y1s, comm); CHKERRQ(ierr);

  ierr = OCE1Set(*oce1, fem, y1s); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode OCE1View(OCE1 self) {
  FEMInfView(self->fem);
  Y1sView(self->y1s);
  return 0;
}

PetscErrorCode OCE1GetSizes(OCE1 self, int *n_r, int *n_y) {
  FEMInfGetSize(self->fem, n_r);
  Y1sGetSize(self->y1s, n_y);
  return 0;
}

PetscErrorCode OCE1SetSr(OCE1 self) {
  FEMInfSetSR1Mat(self->fem, &self->s_r);
  return 0;
}
PetscErrorCode OCE1SetSy(OCE1 self) {
  Y1sSetSY1Mat(self->y1s, &self->s_y);
  return 0;
}

PetscErrorCode OCE1SetSMat(OCE1 self, Mat *M) {

  if(self->s_r == NULL)
    OCE1SetSr(self);
  if(self->s_y == NULL)
    OCE1SetSy(self);

  PetscErrorCode ierr;
  ierr = MatSetSynthesizeFast(self->s_r, self->s_y, self->comm, M);
  CHKERRQ(ierr);

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
  ierr = MatScale(*M, -0.5); CHKERRQ(ierr);
  ierr = MatDestroy(&d2_r); CHKERRQ(ierr);

  Mat l_r, l_y, L;
  ierr = FEMInfSetR2invR1Mat(self->fem, &l_r); CHKERRQ(ierr);
  ierr = Y1sSetLambdaY1Mat(self->y1s, &l_y); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(l_r, l_y, self->comm, &L);
  ierr = MatAXPY(*M, 0.5, L, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&l_r); MatDestroy(&l_y); MatDestroy(&L);

  return 0;
}
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat *M) {

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

      ierr = MatAXPY(*M, zz, V, DIFFERENT_NONZERO_PATTERN);
      
      MatDestroy(&pq_r); MatDestroy(&pq_y); MatDestroy(&V);
    }
  }
  return 0;
  
}
