#include <rescol/mat.h>
#include <rescol/dvr.h>

// ------- Lapack -------
int dgetrf_(long*, long*, PetscReal*,   long*, long*, long*);
int zgetrf_(long*, long*, PetscScalar*, long*, long*, long*);
int dgetri_(long*, double*, long*, long*, double*, long*, long*);
int zgetri_(long*, PetscScalar*, long*, long*, PetscScalar*, long*, long*);

// ------- external functions ---------
PetscErrorCode MatCreateTransformMat(PetscReal *ws, int nq, int ne, 
				     MPI_Comm comm, Mat *A) {
  MatCreate(comm, A);
  MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, nq*ne, (nq-2)*ne + ne-1);
  MatSetUp(*A);

  int ii = 1;
  int jj = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 1; i_q < nq-1; i_q++) {
      MatSetValue(*A, ii, jj, 1.0/sqrt(ws[i_ele*nq+i_q]), INSERT_VALUES);
      ii++; jj++;
    }
    if(i_ele != ne-1) {
      PetscScalar d = 1.0/sqrt(ws[i_ele*nq+nq-1] + ws[(1+i_ele)*nq]);
      MatSetValue(*A, ii, jj, d, INSERT_VALUES);
      ii++;
      MatSetValue(*A, ii, jj, d, INSERT_VALUES);
      ii++; jj++;
    }
  }
  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
  return 0;
}
int NumDVRBasis(int nq, int ne) {
  return ne * (nq-2) + (ne-1);  
}
PetscErrorCode DerivLS(PetscReal *xs, PetscReal *ws, 
		       PetscInt ne, PetscInt nq, 
		       PetscInt i, PetscInt m, PetscInt mp, PetscScalar *v) {
  /*
    gives derivative of Lobatto shape function f_{i,m}
    at quadrature point x_{i,mp}

    nq : number of quadrature in each element
    ne : number of element
    i : index of element
    m,mp : index of quadrature related to Lobatto shape function
   */

  if(!(0 <= i && i < ne)) {
    *v = 0.0;
    return 0;
  }
  
  if(m == mp) {
    if(m == nq-1) 
      *v = +1.0/(2.0*ws[i*nq+m]);
    else if(m == 0) 
      *v = -1.0/(2.0*ws[i*nq+m]);
    else
      *v = 0.0;
  } else {
    *v = 1.0/(xs[i*nq+m]-xs[i*nq+mp]);
    for (int k = 0; k < nq; k++) 
      if(k != m && k != mp) 
	*v *= (xs[i*nq+mp]-xs[i*nq+k])/(xs[i*nq+m]-xs[i*nq+k]);
  }
  return 0;
}

// ------- Basic Methods ---------
PetscErrorCode DVRCreate(MPI_Comm comm, DVR *p_self) {
  DVR self;
  PetscMalloc1(1, &self);
  self->comm = comm;
  self->bps = NULL;
  *p_self = self;
  return 0;
}
PetscErrorCode DVRDestroy(DVR *p_self) {
  PetscErrorCode ierr;
  DVR self = *p_self;
  ierr = BPSDestroy(&self->bps);  CHKERRQ(ierr);    
  ierr = PetscFree(self->xs); CHKERRQ(ierr);    
  ierr = PetscFree(self->ws); CHKERRQ(ierr);    
  ierr = PetscFree(self->xs_basis); CHKERRQ(ierr);    

  if(self->D2_R1LSMat != NULL) {
    ierr = MatDestroy(&self->D2_R1LSMat); CHKERRQ(ierr);    
  }
  if(self->R2_R1LSMat != NULL) {
    ierr = MatDestroy(&self->R2_R1LSMat);CHKERRQ(ierr);  
  }
  if(self->T != NULL) {
    ierr = MatDestroy(&self->T);CHKERRQ(ierr);  
  }
  if(self->TT != NULL) {
    ierr = MatDestroy(&self->TT);CHKERRQ(ierr);  
  }
  if(self->T2 != NULL) {
    ierr = MatDestroy(&self->T2);CHKERRQ(ierr);  
  }
  if(self->T2T != NULL) {
    ierr = MatDestroy(&self->T2T);CHKERRQ(ierr);  
  }
    
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode DVRView(DVR self, PetscViewer v) {

  PetscViewerType type;
  PetscViewerGetType(v, &type);

  if(strcmp(type, "ascii") != 0) 
    SETERRQ(self->comm, 1, "unsupported type");

  PetscViewerASCIIPrintf(v, ">>>> FEM-DVR >>>>\n");
  PetscViewerASCIIPrintf(v, "nq: %d\n", self->nq);  
  PetscViewerASCIIPrintf(v, "num_basis: %d\n", self->num_basis);
  BPSView(self->bps, v);  
  //  ScalerView(self->scaler, v);
  PetscViewerASCIIPrintf(v, "<<<< FEM-DVR <<<<\n");

  return 0;
}

PetscErrorCode DVRSetKnots(DVR self, int nq, BPS bps) {

  PetscErrorCode ierr;

  // data num
  self->nq = nq;
  self->bps = bps;

  int ne; BPSGetNumEle(bps, &ne);
  int nx = nq * ne;
  ierr = PetscMalloc1(nx, &self->xs); CHKERRQ(ierr);
  ierr = PetscMalloc1(nx, &self->ws); CHKERRQ(ierr);
  self->num_basis = NumDVRBasis(self->nq, ne);  
  ierr = PetscMalloc1(self->num_basis, &self->xs_basis); CHKERRQ(ierr);
  
  PetscReal *zs; BPSGetZs(bps, &zs, NULL);
  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    PetscScalar a = zs[i_ele];
    PetscScalar b = zs[i_ele+1];
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w;
      LobGauss(self->nq, i_q, &x, &w);
      x = (b+a)/2.0 + (b-a)/2.0*x;
      w = (b-a)/2.0*w;
      self->xs[idx] = x; self->ws[idx] = w;
      idx++;
    }
  }
  idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 1; i_q < nq-1; i_q++) {
      self->xs_basis[idx] = self->xs[i_ele*nq+i_q];      
      idx++;      
    }
    if(i_ele != ne-1) {
      self->xs_basis[idx] = self->xs[(i_ele+1)*nq];
      idx++;
    }
  }
  self->D2_R1LSMat = NULL;
  self->R2_R1LSMat = NULL;
  ierr = MatCreateTransformMat(self->ws, nq, ne, self->comm, &self->T); 
  CHKERRQ(ierr);
  ierr = MatTranspose(self->T, MAT_INITIAL_MATRIX, &self->TT); CHKERRQ(ierr);
  self->T2 = NULL;
  self->T2T = NULL;

  return 0;
}
PetscErrorCode DVRSetUp(DVR self) {
  return 0;
}
PetscErrorCode DVRSetFromOptions(DVR self) {

  PetscBool find;
  PetscInt nq = 2;  
  PetscErrorCode ierr;
  
  BPS bps;
  BPSCreate(self->comm, &bps);
  ierr = BPSSetFromOptions(self->bps); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-dvr_nq", &nq, &find); CHKERRQ(ierr);
  CHKERRQ(ierr);

  ierr = DVRSetKnots(self, nq, bps); CHKERRQ(ierr);
  return 0;
}

// ---- Accessor ----
PetscErrorCode DVRBasisPsi(DVR self, int i, PetscScalar x, PetscScalar *y) {
  SETERRQ(self->comm, 1, "not implemented yet");
  *y = i + x;
  return 0;
}
PetscErrorCode DVRGetSize(DVR self, int *n) {
  *n = self->num_basis;
  return 0;
}

// ---- inner ----
PetscErrorCode DVRPrepareT2(DVR self) {
  PetscErrorCode ierr;
  ierr = MatSynthesize(self->T, self->T, 1.0, MAT_INITIAL_MATRIX, &self->T2); 
  CHKERRQ(ierr);  
  ierr = MatTranspose(self->T2, MAT_INITIAL_MATRIX, &self->T2T); 
  CHKERRQ(ierr);  
  return 0;
}

// ------- R1Mat/R2Mat ----------
PetscErrorCode DVRCreateR1Mat(DVR self, Mat *M) {
  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int n = (nq-2) * ne + ne -1;
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(*M);
  return 0;
}

PetscErrorCode DVRSR1Mat(DVR self, Mat M) {
  PetscErrorCode ierr;
  
  for(int i = 0; i < self->num_basis; i++) 
    ierr = MatSetValue(M, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRD2R1Mat(DVR self, Mat M) {
  Mat LSR1, T;
  DVRCreateR1LSMat(self, &LSR1);
  //  DVRCreateR1Mat(self, &T);
  DVRD2R1LSMat(self, LSR1);
  DVRLSMatToMat(self, LSR1, MAT_INITIAL_MATRIX, &T);
  MatDestroy(&LSR1);
  MatCopy(T, M, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&T);
  return 0;
}
PetscErrorCode DVRR2invR1Mat(DVR self, Mat M) {
  PetscErrorCode ierr;
  
  for(int i = 0; i < self->num_basis; i++) {
    PetscScalar x = self->xs_basis[i];
    ierr = MatSetValue(M, i, i, 1.0/(x*x), INSERT_VALUES); CHKERRQ(ierr);
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRENR1Mat(DVR self, int q, double a, Mat M) {

  PetscErrorCode ierr;
  
  for(int i = 0; i < self->num_basis; i++) {
    PetscReal x = self->xs_basis[i];
    PetscScalar g = x > a ? x : a;
    PetscScalar s = x < a ? x : a;
    PetscScalar v;
    v = pow(s/g, q)/g;
    ierr = MatSetValue(M, i, i, v, INSERT_VALUES); CHKERRQ(ierr);

    /*
    PetscScalar v;
    PartialCoulomb(q, a, self->xs_basis[i], &v);
    ierr = MatSetValue(*M, i, i, v, INSERT_VALUES); CHKERRQ(ierr);
    */
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVREER2Mat(DVR self, int q, Mat M) {
  PetscErrorCode ierr;
  Mat LS;
  ierr = DVRCreateR1LSMat(self, &LS); CHKERRQ(ierr);
  ierr = DVREER2LSMat(self, q, LS); CHKERRQ(ierr);
  ierr = DVRR2LSMatToR2Mat(self, LS, &M); CHKERRQ(ierr);
  return 0;
}

// ---- LSR1Mat/LSR2Mat ----
PetscErrorCode DVRCreateR1LSMat(DVR self, Mat *M) {
  
  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int n = nq * ne;
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode DVRCreateR2LSMat(DVR self, Mat *M) {

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);

  int n = pow(nq * ne, 2);
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(*M);
  return 0;

}

PetscErrorCode DVRSR1LSMat(DVR self, Mat M) {

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int n = nq * ne;

  for(int k = 0; k < n; k++)
    MatSetValue(M, k, k, self->ws[k], INSERT_VALUES);

  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRPrepareD2R1LSMat(DVR self) {

  Mat M;
  PetscErrorCode ierr;
  PetscScalar *ds_list;
  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  PetscMalloc1(nq*ne*nq, &ds_list);
  //  ds_list = (PetscScalar*)malloc(sizeof(PetscScalar)*nq*ne*nq);
  for(int i = 0; i < ne; i++)
    for(int m = 0; m < nq; m++)
      for(int mp = 0; mp < nq; mp++) {
	int idx = i*nq*nq + m*nq + mp;
	DerivLS(self->xs, self->ws, ne, nq, i, m, mp, &ds_list[idx]);
      }

  DVRCreateR1LSMat(self, &M);

  for(int i = 0; i < ne; i++) 
    for(int ma = 0; ma < nq; ma++)
      for(int mb = 0; mb < nq; mb++) {
	int idx_a = i*nq + ma;
	int idx_b = i*nq + mb;
	PetscScalar v = 0.0;
	for(int k = 0; k < nq; k++) 
	  v -= ds_list[idx_a*nq+k] * ds_list[idx_b*nq+k] * self->ws[i*nq+k];
	MatSetValue(M, idx_a, idx_b, v, INSERT_VALUES);
      }

  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  PetscFree(ds_list);
  self->D2_R1LSMat = M;
  return 0;  
}
PetscErrorCode DVRD2R1LSMat(DVR self, Mat M){
  
  if(self->D2_R1LSMat == NULL) {
    PetscErrorCode ierr;
    ierr = DVRPrepareD2R1LSMat(self); CHKERRQ(ierr);
  }
  
  PetscErrorCode ierr;
  //ierr = MatConvert(self->D2_R1LSMat, MATSAME, MAT_INITIAL_MATRIX, &M); 
  ierr = MatCopy(self->D2_R1LSMat, M, DIFFERENT_NONZERO_PATTERN); 
  CHKERRQ(ierr);  
  return 0;
  
}
PetscErrorCode DVRPrepareR2invR1LSMat(DVR self) {

  Mat M;
  PetscErrorCode ierr;
  DVRCreateR1LSMat(self, &M);

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);

  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) 
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w, v;
      x = self->xs[idx];
      w = self->ws[idx];
      if(i_ele == 0 && i_q == 0)
	v = 0.0;
      else
	v = 1.0/(x*x);

      ierr = MatSetValue(M, idx, idx, v*w, INSERT_VALUES); CHKERRQ(ierr);
      idx++;
    }    
  
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  self->R2_R1LSMat = M;

  return 0;  
}
PetscErrorCode DVRR2invR1LSMat(DVR self, Mat M) {

  if(self->R2_R1LSMat == NULL) {
    PetscErrorCode ierr;
    ierr = DVRPrepareR2invR1LSMat(self); CHKERRQ(ierr);
  }
  
  PetscErrorCode ierr;
  ierr = MatConvert(self->R2_R1LSMat, MATSAME, MAT_INITIAL_MATRIX, &M); 
  CHKERRQ(ierr);  
  return 0;
  
}
PetscErrorCode DVRENR1LSMat(DVR self, int q, PetscReal a, Mat M) {

  PetscErrorCode ierr;
  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);

  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscReal x, w, v;
      x = self->xs[idx];
      w = self->ws[idx];
      if(i_ele == 0 && i_q == 0)
	v = 0.0;
      else {
	PetscScalar g = x > a ? x : a;
	PetscScalar s = x < a ? x : a;
	v = pow(s/g, q)/g;
	// PartialCoulomb(q, a, x, &v);
      }
      ierr = MatSetValue(M, idx, idx, v*w, INSERT_VALUES); CHKERRQ(ierr);
      idx++;
    }    
  }  

  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  
  return 0;
}

PetscErrorCode DVRLSMatToMat(DVR self, Mat A, MatReuse s, Mat *B) {

  MatMatMatMult(self->TT, A, self->T, s, PETSC_DEFAULT, B);
  return 0;  
  
}
PetscErrorCode DVREER2LSMat(DVR self, int q, Mat M) {

  PetscErrorCode ierr;

  Mat D2; DVRCreateR1LSMat(self, &D2);
  Mat R2; DVRCreateR1LSMat(self, &R2);
  ierr = DVRD2R1LSMat(self, D2); CHKERRQ(ierr);
  ierr = DVRR2invR1LSMat(self, R2); CHKERRQ(ierr);
  
  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  PetscReal rmax; BPSGetZMax(self->bps, &rmax);
  
  Mat T;
  ierr = MatConvert(D2, MATSAME, MAT_INITIAL_MATRIX, &T); CHKERRQ(ierr);  
  ierr = MatScale(T, -1.0); CHKERRQ(ierr);  
  if(q != 0)
    ierr = MatAXPY(T, q*(q+1), R2, SUBSET_NONZERO_PATTERN); CHKERRQ(ierr);  
  
  PetscScalar *vs; PetscMalloc1(nq*nq, &vs);
  PetscInt *idx; PetscMalloc1(nq, &idx);

  PetscScalar *work; PetscMalloc1(nq, &work);
  long *ipiv; PetscMalloc1(nq, &ipiv);
  
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int iq = 0; iq < nq; iq++) 
      idx[iq] = i_ele*nq + iq;
    ierr = MatGetValues(T, nq, idx, nq, idx, vs); CHKERRQ(ierr);    
    long info;
    long lnq = nq;
#if defined(PETSC_USE_COMPLEX)
    zgetrf_(&lnq, &lnq, vs, &lnq, ipiv, &info);
    zgetri_(&lnq, vs, &lnq, ipiv, work, &lnq, &info);
#else
    dgetrf_(&lnq, &lnq, vs, &lnq, ipiv, &info);
    dgetri_(&lnq, vs, &lnq, ipiv, work, &lnq, &info);
#endif

    for(int j_ele = 0; j_ele < ne; j_ele++) {
      for(int i = 0; i < nq; i++) 
	for(int j = 0; j < nq; j++) {
	  double ri = self->xs[i_ele*nq+i]; double wi = self->ws[i_ele*nq+i];
	  double rj = self->xs[j_ele*nq+j]; double wj = self->ws[j_ele*nq+j];
	  int idxA = (i_ele*nq+i)*nq*ne + (j_ele*nq+j);
	  double v = pow(ri, q) * pow(rj, q) / pow(rmax, 2*q+1);
	  if(i_ele == j_ele) 
	    if(i != 0 && j != 0)
	      v += (2*q+1) * vs[i*nq+j] / (ri*rj*sqrt(wi)*sqrt(wj));
	  ierr = MatSetValue(M, idxA, idxA, v, INSERT_VALUES); CHKERRQ(ierr);	  
	}
    }
  }

  ierr = MatDestroy(&T); CHKERRQ(ierr);  
  PetscFree(vs); PetscFree(idx); PetscFree(work); PetscFree(ipiv);

  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  
  ierr = MatAssemblyEnd(  M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  
  
  return 0;
}
PetscErrorCode DVRR2LSMatToR2Mat(DVR self, Mat A, Mat *B) {

  PetscErrorCode ierr;

  if(self->T2 == NULL || self->T2T) 
    DVRPrepareT2(self);
  
  ierr = MatMatMatMult(self->T2T, A, self->T2, 
		       MAT_INITIAL_MATRIX, PETSC_DECIDE, B); CHKERRQ(ierr);  
  return 0;
}

