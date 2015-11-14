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
    PetscViewerASCIIPrintf(v, "OCE1 object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "mu: %f\n", self->mu);
    FEMInfView(self->fem, v);
    Y1sView(self->y1s, v);

    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}
PetscErrorCode OCE1ViewFunc(OCE1 self, Vec c, ViewerFunc v) {

  PetscErrorCode ierr;
  PetscViewerType type;
  PetscViewerGetType(v->base, &type);
  if(strcmp(type, "ascii") != 0) {
    char msg[100]; sprintf(msg, "unsupported type: %s", type);
    SETERRQ(self->comm, 1, msg);
  }

  PetscInt num_xs; 
  PetscReal *xs;
  ViewerFuncGetXs(v, &num_xs, &xs);

  int num_r, num_y;
  OCE1GetSizes(self, &num_r, &num_y);


  Vec *cs; ierr = VecGetSplit(c, num_y, &cs); CHKERRQ(ierr);
  for(int i = 0; i < num_xs; i++) {
    PetscReal x = xs[i];
    PetscViewerASCIIPrintf(v->base, "%f ", xs[i]);
    for(int j = 0; j< num_y; j++ ) {
      PetscScalar y;
      FEMInfPsi(self->fem, cs[j], x, &y);
#if defined(PETSC_USE_COMPLEX)
      PetscReal re = PetscRealPart(y);
      PetscReal im = PetscImaginaryPart(y);
      PetscViewerASCIIPrintf(v->base, "%f %f ", re, im);
#else
      PetscViewerASCIIPrintf(v->base, "%f ", y);
#endif
    }
    PetscViewerASCIIPrintf(v->base, "\n");
  }

  ierr = VecRestoreSplit(c, num_y, &cs); CHKERRQ(ierr);

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

PetscErrorCode OCE1CalcSr(OCE1 self) {
  PetscErrorCode ierr;
  ierr = FEMInfCreateMat(self->fem, 1, &self->s_r);CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(self->fem, self->s_r);CHKERRQ(ierr);
  return 0;
}
PetscErrorCode OCE1CalcSy(OCE1 self) {
  PetscErrorCode ierr;
  ierr = Y1sCreateY1Mat(self->y1s, &self->s_y); CHKERRQ(ierr);
  ierr = Y1sSY1Mat(self->y1s, self->s_y);CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode OCE1CreateMat(OCE1 self, Mat *M) {

  PetscErrorCode ierr;
  if(self->s_r == NULL) {
    ierr = OCE1CalcSr(self); CHKERRQ(ierr);
  }
  if(self->s_y == NULL) {
    ierr = OCE1CalcSy(self); CHKERRQ(ierr);
  }

  ierr = MatMatSynthesizeSymbolic(self->s_r, self->s_y, M);
  
  return 0;
}
PetscErrorCode OCE1SMat(OCE1 self, MatReuse scall, Mat *M, PetscBool *is_id) {

  PetscErrorCode ierr;
  PetscBool _is_id;
  FEMInfGetOverlapIsId(self->fem, &_is_id);
  if(is_id != NULL)
    *is_id = _is_id;

  if(_is_id) 
    return 0;

  if(self->s_r == NULL) {
    ierr = OCE1CalcSr(self); CHKERRQ(ierr);
  }
  if(self->s_y == NULL) {
    ierr = OCE1CalcSy(self); CHKERRQ(ierr);
  }

  ierr = MatMatSynthesize(self->s_r, self->s_y, 1.0, scall, M);
  CHKERRQ(ierr);

  return 0;
}
PetscErrorCode OCE1TMat(OCE1 self, MatReuse scall, Mat *M) {

  if(self->s_r == NULL)
    OCE1CalcSr(self);
  if(self->s_y == NULL)
    OCE1CalcSy(self);

  PetscErrorCode ierr;

  Mat d2_r; 
  ierr = FEMInfCreateMat(self->fem, 1, &d2_r); CHKERRQ(ierr);
  ierr = FEMInfD2R1Mat(self->fem, d2_r);CHKERRQ(ierr);
  ierr = MatMatSynthesize(d2_r, self->s_y, -0.5/self->mu, 
			  scall, M); CHKERRQ(ierr);
  ierr = MatDestroy(&d2_r); CHKERRQ(ierr);

  

  Mat l_r;
  ierr = FEMInfCreateMat(self->fem, 1, &l_r); CHKERRQ(ierr);
  ierr = FEMInfR2invR1Mat(self->fem, l_r);    CHKERRQ(ierr);

  Mat l_y;
  ierr = Y1sCreateY1Mat(self->y1s, &l_y); CHKERRQ(ierr);
  ierr = Y1sLambdaY1Mat(self->y1s, l_y); CHKERRQ(ierr);

  Mat L;
  ierr = MatMatSynthesize(l_r, l_y, 0.5/self->mu, 
			  MAT_INITIAL_MATRIX, &L); CHKERRQ(ierr);
  MatDestroy(&l_r); MatDestroy(&l_y); 

  ierr = MatAXPY(*M, 1.0, L, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&L);

  return 0;
}
PetscErrorCode OCE1PotMat(OCE1 self, RotSym sym, Pot pot, MatReuse scall, Mat *M) {

  if(sym != ROT_SCALAR)
    SETERRQ(self->comm, 1, "now only sym=ROT_SCALAR is supported");

  PetscErrorCode ierr;
  if(self->s_y == NULL) {
    ierr = OCE1CalcSy(self); CHKERRQ(ierr);
  }

  Mat r1; 
  ierr = FEMInfCreateMat(self->fem, 1, &r1); CHKERRQ(ierr);
  ierr = FEMInfPotR1Mat(self->fem, pot, r1); CHKERRQ(ierr);

  ierr = MatMatSynthesize(r1, self->s_y, 1.0, scall, M);

  MatDestroy(&r1);
  return 0;
}
PetscErrorCode OCE1PlusPotMat(OCE1 self, RotSym sym, Pot pot, Mat M) {

  PetscErrorCode ierr;

  Mat V;
  ierr = OCE1PotMat(self, sym, pot, MAT_INITIAL_MATRIX, &V); CHKERRQ(ierr);
  ierr = MatAXPY(M, 1.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V);
  return 0;
}
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat M) {

  PetscErrorCode ierr;

  int lmax; Y1sGetMaxL(self->y1s, &lmax);
  int qmax = 2*lmax;
  PetscReal zz = -2.0*z;

  for(int q = 0; q < qmax; q++) {
    PetscBool non0;
    Mat pq_y; 
    ierr = Y1sCreateY1Mat(self->y1s, &pq_y); CHKERRQ(ierr);
    ierr = Y1sPqY1Mat(self->y1s, q, pq_y, &non0); CHKERRQ(ierr);
    if(non0) {
      Mat pq_r;
      ierr = FEMInfCreateMat(self->fem, 1, &pq_r); CHKERRQ(ierr);
      ierr = FEMInfENR1Mat(self->fem, q, a, pq_r); CHKERRQ(ierr);

      Mat V;
      ierr = MatMatSynthesize(pq_r, pq_y, 1.0, 
			      MAT_INITIAL_MATRIX, &V); CHKERRQ(ierr);

      ierr = MatAXPY(M, zz, V, DIFFERENT_NONZERO_PATTERN);
      CHKERRQ(ierr);
      
      
      MatDestroy(&pq_r); MatDestroy(&V);
    }
    MatDestroy(&pq_y);
  }
  return 0;
  
}
