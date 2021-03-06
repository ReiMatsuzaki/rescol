#include <float.h>
#include "../include/math.h"
#include "../include/oce1.h"

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
    PetscViewerASCIIPrintf(v, "fem:");
    FEMInfView(self->fem, v);
    PetscViewerASCIIPrintf(v, "y1s:");
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
  ViewerFuncGetRange(v, &num_xs, &xs);

  int num_r, num_y;
  OCE1GetSizes(self, &num_r, &num_y);

  Vec *cs; ierr = VecGetSplit(c, num_y, &cs); CHKERRQ(ierr);

  PetscViewerASCIIPrintf(v->base, "r,");
  for(int j = 0; j < num_y; j++) {
    int L = self->y1s->ls[j];
    PetscViewerASCIIPrintf(v->base, "re_%d,im_%d", L, L);
    if(j != num_y-1) 
	PetscViewerASCIIPrintf(v->base, ",");
  }
  PetscViewerASCIIPrintf(v->base, "\n");
  
  for(int i = 0; i < num_xs; i++) {
    PetscReal x = xs[i];
    PetscViewerASCIIPrintf(v->base, "%f,", xs[i]);
    for(int j = 0; j< num_y; j++ ) {
      PetscScalar y;
      FEMInfPsiOne(self->fem, cs[j], x, &y);
#if defined(PETSC_USE_COMPLEX)
      PetscReal re = PetscRealPart(y);
      PetscReal im = PetscImaginaryPart(y);
      PetscViewerASCIIPrintf(v->base, "%f,%f", re, im);
#else
      PetscViewerASCIIPrintf(v->base, "%f", y);
#endif
      if(j != num_y-1) 
	PetscViewerASCIIPrintf(v->base, ",");
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

  PetscOptionsGetReal(NULL, NULL, "-oce1_mu", &mu, &find); 
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

PetscErrorCode OCE1Fit(OCE1 self, PF pf, int L, KSP ksp, Vec c) {
  // >> Start LogEvent >>
  int EVENT_id;
  PetscLogEventRegister("OCE1Fit", 0, &EVENT_id);  
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  

  PetscErrorCode ierr;
  if(self->y1s->ls[0] != L) {
    char msg[1000];
    sprintf(msg, "now only fitting for lowest angular momentum is supported.\nL = %d\nlowest L = %d", L, self->y1s->ls[0]);
    SETERRQ(self->comm, 1, msg);
  }

  ierr = VecSet(c, 0.0); CHKERRQ(ierr);
  Vec c_fem;
  ierr = FEMInfCreateVec(self->fem, 1, &c_fem);  CHKERRQ(ierr);
  ierr = FEMInfFit(self->fem, pf, ksp, c_fem);  CHKERRQ(ierr);

  int numr;
  ierr = FEMInfGetSize(self->fem, &numr); CHKERRQ(ierr);
  PetscInt *indices;
  PetscScalar *values;
  ierr = PetscMalloc(numr*(sizeof(PetscInt)), &indices); CHKERRQ(ierr);
  ierr = PetscMalloc(numr*(sizeof(PetscScalar)), &values); CHKERRQ(ierr);
  
  PetscScalar *ptr_c_fem;
  ierr = VecGetArray(c_fem, &ptr_c_fem); CHKERRQ(ierr);
  for(int i = 0; i < numr; i++) {
    indices[i] = i;
    values[i] = ptr_c_fem[i];
  }
  ierr = VecRestoreArray(c_fem, &ptr_c_fem); CHKERRQ(ierr);
  
  ierr = VecSetValues(c, numr, indices, values, INSERT_VALUES); CHKERRQ(ierr);
  
  
  ierr = VecAssemblyBegin(c); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(c); CHKERRQ(ierr);

  ierr = VecDestroy(&c_fem); CHKERRQ(ierr);
  ierr = PetscFree(indices); CHKERRQ(ierr);
  ierr = PetscFree(values); CHKERRQ(ierr);


  // >> End LogEvent >>
  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  
  return 0;
  
}
PetscErrorCode OCE1Psi(OCE1 self, Vec cs, int L, int M, Vec xs, Vec ys) {

  PetscErrorCode ierr;

  if(L != 0 || M != 0)
    SETERRQ(self->comm, 1, "now only L=0 and M=0 is supported");

  int numr;
  ierr = FEMInfGetSize(self->fem, &numr); CHKERRQ(ierr);
    
  PetscInt    *indices;
  PetscScalar *values;
  PetscMalloc(numr*(sizeof(PetscInt)),    &indices);
  PetscMalloc(numr*(sizeof(PetscScalar)), &values);

  PetscScalar *ptr_cs;
  VecGetArray(cs, &ptr_cs);
  for(int i = 0; i < numr; i++) {
    indices[i] = i;
    values[i]  = ptr_cs[i];
  }
  VecRestoreArray(cs, &ptr_cs);

  Vec c_r; FEMInfCreateVec(self->fem, 1, &c_r);
  VecSetValues(c_r, numr, indices, values, INSERT_VALUES);

  FEMInfPsi(self->fem, c_r, xs, ys);
  
  return 0;
}
PetscErrorCode OCE1PsiOne(OCE1 self, Vec cs, int L, int M, PetscScalar x, PetscScalar *y) {

  PetscErrorCode ierr ;

  int numr;
  ierr = FEMInfGetSize(self->fem, &numr); CHKERRQ(ierr);
    
  PetscInt    *indices;
  PetscScalar *values;
  PetscMalloc(numr*(sizeof(PetscInt)),    &indices);
  PetscMalloc(numr*(sizeof(PetscScalar)), &values);

  PetscScalar *ptr_cs;
  ierr =  VecGetArray(cs, &ptr_cs); CHKERRQ(ierr);
  for(int i = 0; i < numr; i++) {
    indices[i] = i;
    values[i]  = ptr_cs[i];
  }
  ierr = VecRestoreArray(cs, &ptr_cs); CHKERRQ(ierr);

  Vec c_r;
  ierr = FEMInfCreateVec(self->fem, 1, &c_r); CHKERRQ(ierr);
  ierr = VecSetValues(c_r, numr, indices, values, INSERT_VALUES); CHKERRQ(ierr);

  ierr = FEMInfPsiOne(self->fem, c_r, x, y); CHKERRQ(ierr);

  PetscFree(indices);
  PetscFree(values);
  VecDestroy(&c_r);
  
  return 0;
}
PetscErrorCode OCE1CalcSr(OCE1 self) {
  int EVENT_id;
  PetscLogEventRegister("OCE1CalcSr", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  
  PetscErrorCode ierr;
  ierr = FEMInfCreateMat(self->fem, 1, &self->s_r);CHKERRQ(ierr);
  ierr = FEMInfSR1Mat(self->fem, self->s_r);CHKERRQ(ierr);

  
  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1CalcSy(OCE1 self) {
  int EVENT_id;
  PetscLogEventRegister("OCE1CalcSy", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  
  PetscErrorCode ierr;
  ierr = Y1sCreateY1Mat(self->y1s, &self->s_y); CHKERRQ(ierr);
  ierr = Y1sSY1Mat(self->y1s, self->s_y);CHKERRQ(ierr);

  
  PetscLogEventEnd(EVENT_id, 0,0,0,0);
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
PetscErrorCode OCE1CreateMatOther(OCE1 self, OCE1 other, Mat *M) {
  
  PetscErrorCode ierr;

  int nr0, ny0, nr1, ny1;
  ierr = OCE1GetSizes(self,  &nr0, &ny0); CHKERRQ(ierr);
  ierr = OCE1GetSizes(other, &nr1, &ny1); CHKERRQ(ierr);

  ierr = MatCreate(self->comm, M);  CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE,
		     nr0*ny0, nr1*ny1); CHKERRQ(ierr);
  ierr = MatSetUp(*M); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode OCE1CreateVec(OCE1 self, Vec *v) {

  int num_r, num_y;
  OCE1GetSizes(self, &num_r, &num_y);
  VecCreate(self->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, num_r*num_y);
  VecSetUp(*v);
  return 0;
  
}
PetscErrorCode OCE1SMat(OCE1 self, MatReuse scall, Mat *M, PetscBool *is_id) {
  int EVENT_id;
  PetscLogEventRegister("OCE1SMat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  

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

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1TMat(OCE1 self, MatReuse scall, Mat *M) {

  int EVENT_id;
  PetscLogEventRegister("OCE1TMat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  
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

  Pot pot_r2;
  ierr = PotCreate(self->comm, &pot_r2); CHKERRQ(ierr);
  ierr = PotSetPower(pot_r2, 1.0, -2); CHKERRQ(ierr);
  Mat l_r;
  ierr = FEMInfCreateMat(self->fem, 1, &l_r); CHKERRQ(ierr);
  //ierr = FEMInfR2invR1Mat(self->fem, l_r);    CHKERRQ(ierr);
  ierr = FEMInfPotR1Mat(self->fem, pot_r2, l_r); CHKERRQ(ierr);

  Mat l_y;
  ierr = Y1sCreateY1Mat(self->y1s, &l_y); CHKERRQ(ierr);
  ierr = Y1sLambdaY1Mat(self->y1s, l_y); CHKERRQ(ierr);

  Mat L;
  ierr = MatMatSynthesize(l_r, l_y, 0.5/self->mu, 
			  MAT_INITIAL_MATRIX, &L); CHKERRQ(ierr);

  ierr = MatAXPY(*M, 1.0, L, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

  MatDestroy(&d2_r);
  PFDestroy(&pot_r2);
  MatDestroy(&l_r);
  MatDestroy(&l_y); 
  MatDestroy(&L);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1PotMat(OCE1 self, RotSym sym, Pot pot, MatReuse scall, Mat *M) {
  int EVENT_id;
  PetscLogEventRegister("OCE1CalcSr", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  

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

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1PlusPotMat(OCE1 self, RotSym sym, Pot pot, Mat M) {
  int EVENT_id;
  PetscLogEventRegister("OCE1CalcSr", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  PetscErrorCode ierr;

  Mat V;
  ierr = OCE1PotMat(self, sym, pot, MAT_INITIAL_MATRIX, &V); CHKERRQ(ierr);
  ierr = MatAXPY(M, 1.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1PlusVneMat(OCE1 self, PetscReal a, PetscReal z, Mat M) {
  int EVENT_id;
  PetscLogEventRegister("OCE1PlusVneMat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  

  PetscErrorCode ierr;

  int lmax; Y1sGetMaxL(self->y1s, &lmax);
  int qmax = 2*lmax + 1;
  PetscReal zz = -2.0*z;

  for(int q = 0; q < qmax; q++) {
    PetscBool non0;
    Mat pq_y; 
    ierr = Y1sCreateY1Mat(self->y1s, &pq_y); CHKERRQ(ierr);
    ierr = Y1sPqY1Mat(self->y1s, q, pq_y, &non0); CHKERRQ(ierr);
    MatAssemblyBegin(pq_y, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pq_y,   MAT_FINAL_ASSEMBLY);
    if(non0) {
      Pot pot;
      ierr = PotCreate(self->comm, &pot); CHKERRQ(ierr);
      ierr = PotSetCoulombNE(pot, q, a, zz); CHKERRQ(ierr);
      
      Mat pq_r;
      ierr = FEMInfCreateMat(self->fem, 1, &pq_r);CHKERRQ(ierr);
      ierr = FEMInfPotR1Mat(self->fem, pot, pq_r); CHKERRQ(ierr);

      Mat V;
      ierr = MatMatSynthesize(pq_r, pq_y, 1.0, 
			      MAT_INITIAL_MATRIX, &V); CHKERRQ(ierr);

      ierr = MatAXPY(M, 1.0, V, DIFFERENT_NONZERO_PATTERN);
      CHKERRQ(ierr);
      
      MatDestroy(&pq_r); MatDestroy(&V);
      PFDestroy(&pot);
    }
    MatDestroy(&pq_y);
  }

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;  
}
PetscErrorCode OCE1ZMat( OCE1 self, OCE1 other, int q, MatReuse scall, Mat *M) {
  int EVENT_id;
  PetscLogEventRegister("OCE1ZMat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  PetscErrorCode ierr;

  if(abs(q) > 1) {
    SETERRQ(self->comm, 1, "q must be -1,0,+1");
  }

  if(self->fem != other->fem) {
    SETERRQ(self->comm, 1, "FEMInf of a and b must be same");
  }

  Pot r1; PotCreate(self->comm, &r1);
  ierr = PotSetPower(r1, 1.0, 1); CHKERRQ(ierr);

  Mat mat_r;
  ierr = FEMInfCreateMat(self->fem, 1, &mat_r); CHKERRQ(ierr);
  ierr = FEMInfPotR1Mat( self->fem, r1, mat_r); CHKERRQ(ierr);

  Mat mat_y;
  ierr = Y1sCreateY1MatOther(self->y1s, other->y1s, &mat_y); CHKERRQ(ierr);
  PetscBool non0;
  ierr = Y1sYqkY1MatOther(self->y1s, other->y1s, 1, q, mat_y, &non0); CHKERRQ(ierr);

  if(!non0) {
    SETERRQ(self->comm, 1, "matrix become empty");
  }
  
  ierr = MatMatSynthesize(mat_r, mat_y, sqrt(4.0*M_PI/3.0), scall, M); CHKERRQ(ierr);
  
  PFDestroy(&r1); MatDestroy(&mat_r); MatDestroy(&mat_y);
  
  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode OCE1DZMat(OCE1 self, OCE1 other, int q, MatReuse scall, Mat *M) {
  int EVENT_id;
  PetscLogEventRegister("OCE1DZMat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  PetscErrorCode ierr;

  if(self->fem != other->fem) {
    SETERRQ(self->comm, 1, "FEMInf of a and b must be same");
  }

  Pot pot_rinv;
  ierr = PotCreate(self->comm, &pot_rinv); CHKERRQ(ierr);
  ierr = PotSetPower(pot_rinv, 1.0, -1); CHKERRQ(ierr);
  
  Mat rinv;
  ierr = FEMInfCreateMat(self->fem, 1, &rinv); CHKERRQ(ierr);
  ierr = FEMInfPotR1Mat(self->fem, pot_rinv, rinv);     CHKERRQ(ierr);
  
  Mat dr;
  ierr = FEMInfCreateMat(self->fem, 1, &dr);  CHKERRQ(ierr);
  ierr = FEMInfD1R1Mat(self->fem, dr);        CHKERRQ(ierr);

  Mat y_dr, y_p, y_m;
  ierr = Y1sCreateY1MatOther(self->y1s, other->y1s, &y_dr); CHKERRQ(ierr);
  ierr = Y1sCreateY1MatOther(self->y1s, other->y1s, &y_p); CHKERRQ(ierr);
  ierr = Y1sCreateY1MatOther(self->y1s, other->y1s, &y_m); CHKERRQ(ierr);
  
  int n_self, n_other;
  ierr = Y1sGetSize(self->y1s, &n_self); CHKERRQ(ierr);
  ierr = Y1sGetSize(other->y1s, &n_other); CHKERRQ(ierr);

  PetscBool find_non0 = PETSC_FALSE;
  for(int i = 0; i < n_self; i++)
    for(int j = 0; j < n_other; j++) {
      int li, mi, lj, mj;
      ierr = Y1sGetLM(self->y1s,  i, &li, &mi); CHKERRQ(ierr);
      ierr = Y1sGetLM(other->y1s, j, &lj, &mj); CHKERRQ(ierr);
      PetscReal pq_ele = Y1EleYqk(li, 1, lj,
				  mi, q, mj) * sqrt(4.0*M_PI/3.0);
      if(fabs(pq_ele) > 10.0*FLT_EPSILON) {
	find_non0 = PETSC_TRUE;
      } else {
	continue;
      }
      int idxm[1]; idxm[0] = i;
      int idxn[1]; idxn[0] = j;
      PetscScalar vs[1]; 
      if(li==lj+1) {
	vs[0] = pq_ele;
	ierr = MatSetValues(y_dr, 1, idxm, 1, idxn, vs, INSERT_VALUES); CHKERRQ(ierr);
	vs[0] = -1.0*(lj+1)*pq_ele;
	ierr = MatSetValues(y_p, 1, idxm, 1, idxn, vs, INSERT_VALUES); CHKERRQ(ierr);
      } else if(li==lj-1) {
	vs[0] = pq_ele;
	ierr = MatSetValues(y_dr, 1, idxm, 1, idxn, vs, INSERT_VALUES); CHKERRQ(ierr);
	vs[0] = 1.0*lj*pq_ele;
	ierr = MatSetValues(y_m, 1, idxm, 1, idxn, vs, INSERT_VALUES); CHKERRQ(ierr);
      }
    }

  if(!find_non0) {
    SETERRQ(self->comm, 1, "matrix becomes empty");
  }

  MatAssemblyBegin(y_dr,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(y_dr, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(y_p, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(y_p,  MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(y_m, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(y_m,  MAT_FINAL_ASSEMBLY);
  
  Mat B, C;
  ierr = MatMatSynthesize(dr,   y_dr, 1.0, scall, M); CHKERRQ(ierr);
  ierr = MatMatSynthesize(rinv, y_p,  1.0, MAT_INITIAL_MATRIX, &B); CHKERRQ(ierr);
  ierr = MatMatSynthesize(rinv, y_m,  1.0, MAT_INITIAL_MATRIX, &C); CHKERRQ(ierr);
  ierr = MatAXPY(*M, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatAXPY(*M, 1.0, C, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  
  PFDestroy(&pot_rinv); MatDestroy(&rinv); MatDestroy(&dr);  
  MatDestroy(&y_dr); MatDestroy(&y_p); MatDestroy(&y_m);
  MatDestroy(&B); MatDestroy(&C);
  
  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;  
}

// -- not used now
PetscErrorCode OceH2plusMatMultH(Mat H, Vec x, Vec y) {
  PetscErrorCode ierr;
  OceH2plus ctx;
  
  ierr = MatShellGetContext(H, &ctx); CHKERRQ(ierr);

  ierr = MatMatDecomposedMult(ctx->d2_r1, ctx->s_y1, -0.5, x, y); CHKERRQ(ierr);
  ierr = MatMatDecomposedMultAdd(ctx->r2inv_r1, ctx->lambda_y1, 0.5, x, y); CHKERRQ(ierr);

  for(int i = 0; i < ctx->nq; i++) {
    ierr = MatMatDecomposedMultAdd(ctx->ne_r1[i], ctx->pq_y1[i], -2.0, x, y); 
    CHKERRQ(ierr);
  }

  return 0;
}
// -- not used now
PetscErrorCode OceH2plusMatMultS(Mat S, Vec x, Vec y) {
  PetscErrorCode ierr;
  OceH2plus ctx;
  
  ierr = MatShellGetContext(S, &ctx); CHKERRQ(ierr);
  ierr = MatMatDecomposedMult(ctx->s_r1, ctx->s_y1, 1.0, x, y); CHKERRQ(ierr);

  return 0;
}
// -- not used now
PetscErrorCode OCE1CreateH2plus(OCE1 self, PetscReal a, PetscReal z, OceH2plus *p_ctx) {
  PetscErrorCode ierr;
  OceH2plus ctx;

  SETERRQ(self->comm, 1, "not supported now");
  
  ierr = PetscNew(&ctx); CHKERRQ(ierr);

  ctx->a = a; ctx->z = z;

  FEMInfCreateMat(self->fem,1 , &ctx->s_r1); 
  FEMInfCreateMat(self->fem,1 , &ctx->d2_r1); 
  FEMInfCreateMat(self->fem,1 , &ctx->r2inv_r1); 
  Y1sCreateY1Mat(self->y1s , &ctx->s_y1);
  Y1sCreateY1Mat(self->y1s , &ctx->lambda_y1); 

  if(self->s_r == NULL) {
    ierr = OCE1CalcSr(self); CHKERRQ(ierr);
  }
  if(self->s_y == NULL) {
    ierr = OCE1CalcSy(self); CHKERRQ(ierr);
  }

  // calculate ctx 
  MatCopy(self->s_r, ctx->s_r1, DIFFERENT_NONZERO_PATTERN);
  FEMInfD2R1Mat(self->fem, ctx->d2_r1);
  //FEMInfR2invR1Mat(self->fem, ctx->r2inv_r1);  
  MatCopy(self->s_y, ctx->s_y1, DIFFERENT_NONZERO_PATTERN);
  Y1sLambdaY1Mat(self->y1s, ctx->lambda_y1);

  int lmax; Y1sGetMaxL(self->y1s, &lmax);
  int qmax = 2*lmax+1; 
  ierr = PetscMalloc1(qmax, &ctx->q); CHKERRQ(ierr);
  ierr = PetscMalloc1(qmax, &ctx->ne_r1); CHKERRQ(ierr);
  ierr = PetscMalloc1(qmax, &ctx->pq_y1); CHKERRQ(ierr);
  int idx = 0;
  for(int q = 0; q < qmax; q++) {
    Mat pq_y, pq_r;
    PetscBool non0;
    ierr = Y1sCreateY1Mat(self->y1s, &pq_y);      CHKERRQ(ierr);
    ierr = FEMInfCreateMat(self->fem, 1, &pq_r); CHKERRQ(ierr);
    ierr = Y1sPqY1Mat(self->y1s, q, pq_y, &non0); CHKERRQ(ierr);
    if(non0) {
      //ierr = FEMInfENR1Mat(self->fem, q, a, pq_r); CHKERRQ(ierr);
      ctx->q[idx] = q;
      ctx->ne_r1[idx] = pq_r;
      ctx->pq_y1[idx] = pq_y;
      idx++;
    } else {
      MatDestroy(&pq_y);
      MatDestroy(&pq_r);
    }
  }
  ctx->nq = idx;

  *p_ctx = ctx;
  return 0;
}
// -- not used now
PetscErrorCode OCE1H2plusMat(OCE1 self, OceH2plus ctx, Mat *H, Mat *S, PetscBool *is_id) {
			     
  PetscErrorCode ierr;

  // set overlap property
  PetscBool _is_id;
  FEMInfGetOverlapIsId(self->fem, &_is_id);
  if(is_id != NULL)
    *is_id = _is_id;

  // Mat creation
  Mat _H, _S;
  int nr, ny; FEMInfGetSize(self->fem, &nr); Y1sGetSize(self->y1s, &ny);
  int n = nr*ny;
  ierr = MatCreateShell(self->comm, PETSC_DECIDE, PETSC_DECIDE, 
			n, n, ctx, &_H); CHKERRQ(ierr);
  ierr = MatCreateShell(self->comm, PETSC_DECIDE, PETSC_DECIDE, 
			n, n, ctx, &_S);  CHKERRQ(ierr);

  // set context
  ierr = MatShellSetOperation(_H, MATOP_MULT, 
			      (void(*)(void))OceH2plusMatMultH); CHKERRQ(ierr);
  ierr = MatShellSetOperation(_S, MATOP_MULT, 
			      (void(*)(void))OceH2plusMatMultS);CHKERRQ(ierr);
  *H = _H;
  *S = _S;

  return 0;
}
// -- not used now			       
PetscErrorCode OCE1H2plusMat_direct(OCE1 self, OceH2plus ctx, Mat *H, Mat *S, PetscBool *is_id) {
  
  PetscErrorCode ierr;

  PetscReal a = ctx->a;
  PetscReal z = ctx->z;

  Mat _H;
  ierr = OCE1TMat(self, MAT_INITIAL_MATRIX, &_H); CHKERRQ(ierr);
  ierr = OCE1PlusVneMat(self, a, z, _H); CHKERRQ(ierr);
  *H = _H;

  Mat _S;
  ierr = OCE1SMat(self, MAT_INITIAL_MATRIX, &_S, is_id); CHKERRQ(ierr);
  *S = _S;
  return 0;
}
// -- not used now
PetscErrorCode OCE1H2plusDestroy(OceH2plus *p_ctx) {
  PetscErrorCode ierr;
  OceH2plus ctx = *p_ctx;
  ierr = MatDestroy(&ctx->s_r1); CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->d2_r1); CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->r2inv_r1); CHKERRQ(ierr);

  ierr = MatDestroy(&ctx->s_y1); CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->lambda_y1); CHKERRQ(ierr);

  for(int q = 0; q < ctx->nq; q++) {
    ierr = MatDestroy(&ctx->ne_r1[q]); CHKERRQ(ierr);
    ierr = MatDestroy(&ctx->pq_y1[q]); CHKERRQ(ierr);
  }
  ierr = PetscFree(ctx->q); CHKERRQ(ierr);
  ierr = PetscFree(ctx->ne_r1); CHKERRQ(ierr);
  ierr = PetscFree(ctx->pq_y1); CHKERRQ(ierr);

  PetscFree(*p_ctx);

  return 0;
}

/*
PetscErrorCode OceH2plusCreate(OceH2plus *p_self) {
  PetscErrorCode ierr;
  OceH2plus self;  
  ierr = PetscNew(&self); CHKERRQ(ierr);
  *p_self = self;
  return 0;	       
}
PetscErrorCode OceH2plusDestroy(OceH2plus *p_self) {
  PetscErrorCode ierr;
  OceH2plus self = *p_self;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode OceH2plusCalc(OceH2plus self, FEMInf fem, Y1s y1s) {

  FEMInfGetSize(fem, &self->nr);
  Y1sGetSize(y1s, &self->ny);
  FEMInfGetOverlapIsId(fem, &self->is_id);

  FEMInfCreateMat(fem, 1, &self->s_r1);     FEMInfSR1Mat(fem, self->s_r1);
  FEMInfCreateMat(fem, 1, &self->d2_r1);    FEMInfD2R1Mat(fem, self->d2_r1);  
  FEMInfCreateMat(fem, 1, &self->r2inv_r1); FEMInfD2R1Mat(fem, self->r2inv_r1);  

  Y1sCreateY1Mat(y1s, &self->s_y1); Y1sSY1Mat(y1s, self->s_y1);
  Y1sCreateY1Mat(y1s, &self->lambda_y1); Y1sLambdaY1Mat(y1s, self->lambda_y1);  

  return 0;
}
PetscErrorCode OceH2PlusCreateMat(OceH2plus self, Mat *M);
PetscErrorCode OceH2plusHMat(OceH2plus self, Mat H);
PetscErrorCode OceH2plusSMat(OceH2plus self, Mat S, PetscBool *is_id);

PetscErrorCode OceH2plusView(OceH2plus self, PetscViewer v) {
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  
  ierr = OceH2plusCheckState(self); CHKERRQ(ierr);
  
  if(v == NULL)
    PetscViewerASCIIGetStdout(self->comm, &v);
  
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);  
  
  if(iascii) {
    PetscViewerASCIIPrintf(v, "OceH2plus object:\n");
    PetscViewerASCIIPushTab(v);
    
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {
    
  } else if(isdraw) {
    
  }
  return 0;
  
}

*/



