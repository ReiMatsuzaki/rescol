#include <rescol/oce2.h>

PetscErrorCode OCE2Create(MPI_Comm comm, OCE2 *p_self) {

  OCE2 self;
  PetscNew(&self);
  *p_self = NULL;

  self->comm = comm;

  self->fem = NULL;
  self->y2s = NULL;

  self->s_r1 = NULL;
  self->s_y2 = NULL;
  
  *p_self = self;
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
PetscErrorCode OCE2SetFromOptions(OCE2 self) {

  PetscErrorCode ierr;
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetFromOptions(fem);
  Y2s y2s;    Y2sCreate(comm, &y2s);    Y2sSetFromOptions(y2s);

  ierr = OCE2Set(self, fem, y2s); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode OCE2View(OCE2 self, PetscViewer v) {

  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  PetscViewerType type;     PetscViewerGetType(v, &type);
  PetscViewerFormat format; PetscViewerGetFormat(v, &format);

  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);    

  if(iascii) {
    PetscViewerASCIIPrintf(v, "OCE2 object:\n");
    PetscViewerASCIIPushTab(v);
    FEMInfView(self->fem, v);
    Y2sView(self->y2s, v);
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}

PetscErrorCode OCE2GetSizes(OCE2 this, int *n_r1, int *n_y2) {
  FEMInfGetSize(this->fem, n_r1);
  Y2sGetSize(this->y2s, n_y2);
  return 0;
}

PetscErrorCode OCE2SetSr1(OCE2 self) {

  PetscErrorCode ierr;
  ierr = FEMInfCreateMat(self->fem, 1, &self->s_r1); CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(self->fem, self->s_r1); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode OCE2SetSy2(OCE2 self) {
  
  PetscErrorCode ierr;
  ierr = Y2sCreateY2Mat(self->y2s, &self->s_y2); CHKERRQ(ierr);
  ierr = Y2sSY2Mat(self->y2s, self->s_y2); CHKERRQ(ierr);
  CHKERRQ(ierr);
  return 0;

}

PetscErrorCode OCE2CreateMat(OCE2 self, Mat *M) {

  int nr, ny;
  FEMInfGetSize(self->fem, &nr);
  Y2sGetSize(self->y2s, &ny);
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nr*ny, nr*ny);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode OCE2SMat(OCE2 self, MatReuse scall, Mat *M, PetscBool *is_id) {
  
  PetscBool _is_id; 
  FEMInfGetOverlapIsId(self->fem, &_is_id);
  if(is_id != NULL)
    *is_id = _is_id;
  
  if(_is_id) 
    return 0;

  PetscErrorCode ierr;
  if(self->s_r1 == NULL) {
    ierr = OCE2SetSr1(self); CHKERRQ(ierr);
  }
  if(self->s_y2 == NULL) {
    ierr = OCE2SetSy2(self); CHKERRQ(ierr);
  }

  ierr = MatMatMatSynthesize(self->s_r1, self->s_r1, self->s_y2, 1.0, scall, M);
  CHKERRQ(ierr);
  return 0;
 
}
PetscErrorCode OCE2TMat(OCE2 self, MatReuse scall, Mat *M) {

  PetscErrorCode ierr;
  if(self->s_r1 == NULL) {
    ierr = OCE2SetSr1(self); CHKERRQ(ierr);
  }
  if(self->s_y2 == NULL) {
    ierr = OCE2SetSy2(self); CHKERRQ(ierr);
  }

  Mat d2_r;
  ierr = FEMInfCreateMat(self->fem, 1, &d2_r); CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(self->fem, d2_r); CHKERRQ(ierr);
  ierr = MatMatMatSynthesize(d2_r, self->s_r1, self->s_y2, -0.5, scall, M);CHKERRQ(ierr);
  
  Mat D;
  ierr = MatMatMatSynthesize(self->s_r1, d2_r, self->s_y2, -0.5, 
			     MAT_INITIAL_MATRIX, &D);CHKERRQ(ierr);
  ierr = MatAXPY(*M, 1.0, D, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&d2_r); MatDestroy(&D);
  
  Mat l_r1; 
  ierr = FEMInfCreateMat(self->fem, 1, &l_r1); FEMInfR2invR1Mat(self->fem, l_r1);

  Mat l1_y2;
  PetscBool non0;
  ierr = Y2sCreateY2Mat(self->y2s, &l1_y2); Y2sLambda1Y2Mat(self->y2s, l1_y2, &non0);
  if(non0) {
    Mat L;
    ierr = MatMatMatSynthesize(l_r1, self->s_r1, l1_y2, 0.5, MAT_INITIAL_MATRIX, &L);
    ierr = MatAXPY(*M, 1.0, L, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&L);
  }
  MatDestroy(&l1_y2);

  Mat l2_y2;
  ierr = Y2sCreateY2Mat(self->y2s, &l2_y2); Y2sLambda2Y2Mat(self->y2s, l2_y2, &non0);
  if(non0) {
    Mat L;
    ierr = MatMatMatSynthesize(self->s_r1, l_r1, l2_y2, 0.5, MAT_INITIAL_MATRIX, &L);
    ierr = MatAXPY(*M, 1.0, L, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&L);
  }
  MatDestroy(&l1_y2);  
  MatDestroy(&l_r1);
  return 0;
}
PetscErrorCode OCE2PlusVneMat(OCE2 self, PetscReal a, PetscReal z, Mat M) {

  if(self->s_r1 == NULL)
    OCE2SetSr1(self);

  PetscErrorCode ierr;

  int lmax; Y2sGetMaxL(self->y2s, &lmax);
  int qmax = 2*lmax;

  PetscScalar zz = -2.0*z;
  for(int q = 0; q <= qmax; q++) {
    PetscBool non0_1, non0_2;
    Mat pq1A_y2;    
    ierr = Y2sCreateY2Mat(self->y2s, &pq1A_y2); CHKERRQ(ierr);
    ierr = Y2sPq1AY2Mat(self->y2s, q, pq1A_y2, &non0_1); CHKERRQ(ierr);
    Mat pq2A_y2;
    ierr = Y2sCreateY2Mat(self->y2s, &pq2A_y2); CHKERRQ(ierr);
    ierr = Y2sPq1AY2Mat(self->y2s, q, pq2A_y2, &non0_2); CHKERRQ(ierr);
    if(non0_1 && non0_2) {
      Mat q_r1;
      ierr = FEMInfCreateMat(self->fem, 1, &q_r1); CHKERRQ(ierr);
      ierr = FEMInfENR1Mat(self->fem, q, a, q_r1); CHKERRQ(ierr);

      Mat V1;
      ierr = MatMatMatSynthesize(q_r1, self->s_r1, pq1A_y2, 1.0,
				 MAT_INITIAL_MATRIX, &V1); CHKERRQ(ierr);
      ierr = MatAXPY(M, zz, V1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&V1);

      Mat V2;
      ierr = MatMatMatSynthesize(self->s_r1, q_r1, pq2A_y2, 1.0,
				    MAT_INITIAL_MATRIX, &V2); CHKERRQ(ierr);
      ierr = MatAXPY(M, zz, V2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&V2);
      
      MatDestroy(&q_r1);
    }
    MatDestroy(&pq1A_y2); MatDestroy(&pq2A_y2);
  }

  return 0;
}
PetscErrorCode OCE2PlusVeeMat(OCE2 self, Mat M) {
  
  PetscErrorCode ierr;
  
  int lmax; Y2sGetMaxL(self->y2s, &lmax);
  int qmax = 2*lmax;

  for(int q = 0; q < qmax; q++) {
    PetscBool non0;
    Mat pq;
    Y2sCreateY2Mat(self->y2s, &pq);
    Y2sPq12Y2Mat(self->y2s, q, pq, &non0);
    if(non0) {
      Mat r2, V;
      ierr = FEMInfCreateMat(self->fem, 2, &r2); CHKERRQ(ierr);
      ierr = FEMInfEER2Mat(self->fem, q, r2); CHKERRQ(ierr);
      ierr = MatMatSynthesize(r2, pq, 1.0, MAT_INITIAL_MATRIX, &V); CHKERRQ(ierr);
      ierr = MatAXPY(M, 1.0, V, DIFFERENT_NONZERO_PATTERN);
      MatDestroy(&V); MatDestroy(&r2);
    }
    MatDestroy(&pq); 
  }
  return 0;
}



