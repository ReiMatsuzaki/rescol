#include "oce2.h"

PetscErrorCode OCE2Create(OCE2 *oce2, MPI_Comm comm) {

  OCE2 _oce2;
  PetscNew(&_oce2);
  *oce2 = NULL;

  _oce2->comm = comm;

  _oce2->fem = NULL;
  _oce2->y2s = NULL;

  _oce2->s_r1 = NULL;
  _oce2->s_y2 = NULL;
  
  *oce2 = _oce2;
  return 0;
}
PetscErrorCode OCE2Destroy(OCE2 *oce2) {

  PetscErrorCode ierr;
  
  if((*oce2)->fem) {
    ierr = FEMInfDestroy(&(*oce2)->fem); CHKERRQ(ierr);
  }

  if((*oce2)->y2s) {
    ierr = Y2sDestroy(&(*oce2)->y2s); CHKERRQ(ierr);
  }

  if((*oce2)->s_r1) {
    ierr = MatDestroy(&(*oce2)->s_r1); CHKERRQ(ierr);
  }

  if((*oce2)->s_y2) {
    ierr = MatDestroy(&(*oce2)->s_y2); CHKERRQ(ierr);
  }

  ierr = PetscFree(*oce2); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode OCE2Set(OCE2 this, FEMInf fem, Y2s y2s) {

  this->fem = fem;
  this->y2s = y2s;
  return 0;
}
PetscErrorCode OCE2CreateFromOptions(OCE2 *oce2, MPI_Comm comm) {

  PetscErrorCode ierr;
  ierr = OCE2Create(oce2, comm); CHKERRQ(ierr);
  
  FEMInf fem; 
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);

  Y2s y2s;
  ierr = Y2sCreateFromOptions(&y2s, comm); CHKERRQ(ierr);

  ierr = OCE2Set(*oce2, fem, y2s); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode OCE2View(OCE2 this) {

  PetscErrorCode ierr;
  ierr = FEMInfView(this->fem); CHKERRQ(ierr);
  ierr = Y2sView(this->y2s); CHKERRQ(ierr);
  return 0;

}

PetscErrorCode OCE2GetSizes(OCE2 this, int *n_r1, int *n_y2) {
  FEMInfGetSize(this->fem, n_r1);
  Y2sGetSize(this->y2s, n_y2);
  return 0;
}

PetscErrorCode OCE2SetSr1(OCE2 this) {

  PetscErrorCode ierr;
  ierr = FEMInfSetSR1Mat(this->fem, &this->s_r1);
  CHKERRQ(ierr);
  return 0;

}
PetscErrorCode OCE2SetSy2(OCE2 this) {
  
  PetscErrorCode ierr;
  ierr = Y2sSetSY2Mat(this->y2s, &this->s_y2);
  CHKERRQ(ierr);
  return 0;

}

PetscErrorCode OCE2SetSMat(OCE2 this, Mat *M) {

  if(this->s_r1 == NULL)
    OCE2SetSr1(this);

  if(this->s_y2 == NULL)
    OCE2SetSy2(this);

  PetscErrorCode ierr;
  ierr = MatSetSynthesize3Fast(this->s_r1, this->s_r1, this->s_y2,
			       this->comm, M); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode OCE2SetTMat(OCE2 this, Mat *M) {

  if(this->s_r1 == NULL)
    OCE2SetSr1(this);

  if(this->s_y2 == NULL)
    OCE2SetSy2(this);

  PetscErrorCode ierr;

  Mat d2_r1;
  ierr = FEMInfSetD2R1Mat(this->fem, &d2_r1); CHKERRQ(ierr);
  ierr = MatScale(d2_r1, -0.5); CHKERRQ(ierr);
  ierr = MatSetSynthesize3Fast(d2_r1, this->s_r1, this->s_y2, 
			       this->comm, M); CHKERRQ(ierr);

  Mat D;
  ierr = MatSetSynthesize3Fast(this->s_r1, d2_r1, this->s_y2,
			       this->comm, &D); CHKERRQ(ierr);
  ierr = MatAXPY(*M, 1.0, D, DIFFERENT_NONZERO_PATTERN);

  MatDestroy(&D); MatDestroy(&d2_r1);

  
  /*
  Mat l_r1;
  ierr = FEMInfSetR2invR1Mat(this->fem, &l_r1); CHKERRQ(ierr);
  
  Mat l1_y2;
  ierr = Y2sSetLambda1Y2Mat(this->y2s, &l1_y2); CHKERRQ(ierr);
  if(l1_y2) {
    Mat  L1;
    ierr = MatSetSynthesize3Fast(l_r1, this->s_r1, l1_y2, 
				 this->comm, &L1);
    ierr = MatAXPY(*M, 0.5, L1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    MatDestroy(&l1_y2); MatDestroy(&L1);
  }

  Mat l2_y2;
  ierr = Y2sSetLambda2Y2Mat(this->y2s, &l2_y2); CHKERRQ(ierr);
  if(l2_y2) {
    Mat  L2;
    ierr = MatSetSynthesize3Fast(this->s_r1, l_r1, l2_y2, 
				 this->comm, &L2); CHKERRQ(ierr);
    ierr = MatAXPY(*M, 0.5, L2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    MatDestroy(&l2_y2); MatDestroy(&L2);
  }

  MatDestroy(&l_r1);
  */
  return 0;

}
PetscErrorCode OCE2PlusVneMat(OCE2 this, PetscReal a, PetscReal z, Mat *M) {

  if(this->s_r1 == NULL)
    OCE2SetSr1(this);

  PetscErrorCode ierr;

  int lmax; Y2sGetMaxL(this->y2s, &lmax);
  int qmax = 2*lmax;

  PetscScalar zz = -2.0*z;
  for(int q = 0; q <= qmax; q++) {
    Mat pq1A_y2, pq2A_y2;
    ierr = Y2sSetPq1AY2Mat(this->y2s, q, &pq1A_y2); CHKERRQ(ierr);
    ierr = Y2sSetPq2AY2Mat(this->y2s, q, &pq2A_y2); CHKERRQ(ierr);
    
    if(pq1A_y2 != NULL && pq2A_y2 != NULL) {
      Mat q_r1;
      ierr = FEMInfSetENR1Mat(this->fem, q, a, &q_r1); CHKERRQ(ierr);

      Mat V1;
      ierr = MatSetSynthesize3Fast(q_r1, this->s_r1, pq1A_y2,
				   this->comm, &V1); CHKERRQ(ierr);
      ierr = MatAXPY(*M, zz, V1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&V1);

      Mat V2;
      ierr = MatSetSynthesize3Fast(this->s_r1, q_r1, pq2A_y2,
				   this->comm, &V2); CHKERRQ(ierr);
      ierr = MatAXPY(*M, zz, V2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&V2);
      
      MatDestroy(&q_r1); MatDestroy(&pq1A_y2); MatDestroy(&pq2A_y2);
    }
  }

  return 0;
}
PetscErrorCode OCE2PlusVeeMat(OCE2 this, Mat *M) {
  
  PetscErrorCode ierr;
  
  int lmax; Y2sGetMaxL(this->y2s, &lmax);
  int qmax = 2*lmax;

  for(int q = 0; q < qmax; q++) {
    Mat pq;
    Y2sSetPq12Y2Mat(this->y2s, q, &pq);
    if(pq != NULL) {
      Mat r2, V;
      ierr = FEMInfSetEER2Mat(this->fem, q, &r2); CHKERRQ(ierr);
      ierr = MatSetSynthesizeFast(r2, pq, this->comm, &V); CHKERRQ(ierr);
      ierr = MatAXPY(*M, 1.0, V, DIFFERENT_NONZERO_PATTERN);
      MatDestroy(&pq); MatDestroy(&V); MatDestroy(&r2);
    }
  }

  return 0;
}

