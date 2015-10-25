#include "mat.h"
#include "dvr.h"

// ------- Lapack -------
int dgetrf_(long*, long*, double*, long*, long*, long*);
int dgetri_(long*, double*, long*, long*, double*, long*, long*);

// ------- external functions ---------
PetscErrorCode MatCreateTransformMat(PetscScalar *ws, int nq, int ne, 
				     MPI_Comm comm, Mat *A) {
  MatCreate(comm, A);
  MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, nq*ne, (nq-2)*ne + ne-1);
  MatSetFromOptions(*A);
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
PetscErrorCode DerivLS(PetscScalar *xs, PetscScalar *ws, 
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
PetscErrorCode DVRCreate(DVR *dvr, int nq, BPS bps, MPI_Comm comm) {

  PetscErrorCode ierr;
  DVR _dvr;
  PetscMalloc1(1, &_dvr);
  *dvr = NULL;

  // data num
  _dvr->comm = comm;  
  _dvr->nq = nq;
  _dvr->bps = bps;

  int ne; BPSGetNumEle(bps, &ne);
  int nx = nq * ne;
  ierr = PetscMalloc1(nx, &_dvr->xs); CHKERRQ(ierr);
  ierr = PetscMalloc1(nx, &_dvr->ws); CHKERRQ(ierr);
  _dvr->num_basis = NumDVRBasis(_dvr->nq, ne);  
  ierr = PetscMalloc1(_dvr->num_basis, &_dvr->xs_basis); CHKERRQ(ierr);
  
  PetscScalar *zs; BPSGetZs(bps, &zs, NULL);
  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    PetscScalar a = zs[i_ele];
    PetscScalar b = zs[i_ele+1];
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w;
      LobGauss(_dvr->nq, i_q, &x, &w);
      x = (b+a)/2.0 + (b-a)/2.0*x;
      w = (b-a)/2.0*w;
      _dvr->xs[idx] = x; _dvr->ws[idx] = w;
      idx++;
    }
  }
  idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 1; i_q < nq-1; i_q++) {
      _dvr->xs_basis[idx] = _dvr->xs[i_ele*nq+i_q];      
      idx++;      
    }
    if(i_ele != ne-1) {
      _dvr->xs_basis[idx] = _dvr->xs[(i_ele+1)*nq];
      idx++;
    }
  }
  _dvr->D2_R1LSMat = NULL;
  _dvr->R2_R1LSMat = NULL;
  ierr = MatCreateTransformMat(_dvr->ws, nq, ne, comm, &_dvr->T); 
  CHKERRQ(ierr);
  ierr = MatTranspose(_dvr->T, MAT_INITIAL_MATRIX, &_dvr->TT); CHKERRQ(ierr);
  _dvr->T2 = NULL;
  _dvr->T2T = NULL;

  *dvr = _dvr;
  return 0;
}
PetscErrorCode DVRCreateFromOptions(DVR *dvr, MPI_Comm comm) {
  PetscBool find;
  BPS bps;
  PetscInt nq = 2;  
  PetscErrorCode ierr;

  ierr = BPSCreate(&bps, comm); CHKERRQ(ierr);
  ierr = BPSSetFromOptions(bps); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-dvr_nq", &nq, &find); CHKERRQ(ierr);
  CHKERRQ(ierr);

  ierr = DVRCreate(dvr, nq, bps, comm); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode DVRDestroy(DVR *dvr) {
  PetscErrorCode ierr;
  DVR this = *dvr;
  ierr = BPSDestroy(&this->bps);  CHKERRQ(ierr);    
  ierr = PetscFree(this->xs); CHKERRQ(ierr);    
  ierr = PetscFree(this->ws); CHKERRQ(ierr);    
  ierr = PetscFree(this->xs_basis); CHKERRQ(ierr);    

  if(this->D2_R1LSMat != NULL) {
    ierr = MatDestroy(&this->D2_R1LSMat); CHKERRQ(ierr);    
  }
  if(this->R2_R1LSMat != NULL) {
    ierr = MatDestroy(&this->R2_R1LSMat);CHKERRQ(ierr);  
  }
  if(this->T != NULL) {
    ierr = MatDestroy(&this->T);CHKERRQ(ierr);  
  }
  if(this->TT != NULL) {
    ierr = MatDestroy(&this->TT);CHKERRQ(ierr);  
  }
  if(this->T2 != NULL) {
    ierr = MatDestroy(&this->T2);CHKERRQ(ierr);  
  }
  if(this->T2T != NULL) {
    ierr = MatDestroy(&this->T2T);CHKERRQ(ierr);  
  }
    
  ierr = PetscFree(*dvr); CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode DVRFPrintf(DVR this, FILE *file, int lvl) {
  
  PetscErrorCode ierr;
  
  if(lvl != 0) {
    SETERRQ(this->comm, 1, "now only lvl=0 is supported");
  }

  PetscFPrintf(this->comm, file, "===== Begin DVR ======\n");
  ierr = BPSFPrintf(this->bps, file, lvl); CHKERRQ(ierr);
  PetscFPrintf(this->comm, file, "#of quad: %d\n", this->nq);
  PetscFPrintf(this->comm, file, "===== End DVR ======\n");
  return 0;
}
PetscErrorCode DVRBasisPsi(DVR this, int i, PetscScalar x, PetscScalar *y) {
  SETERRQ(this->comm, 1, "not implemented yet");
  return 0;
}

// ------- inner ----------
PetscErrorCode DVRPrepareT2(DVR this) {
  PetscErrorCode ierr;
  ierr = MatSetSynthesizeFast(this->T, this->T, this->comm, &this->T2); 
  CHKERRQ(ierr);  
  ierr = MatTranspose(this->T2, MAT_INITIAL_MATRIX, &this->T2T); 
  CHKERRQ(ierr);  
  return 0;
}

// ------- R1Mat ----------
PetscErrorCode DVRInitR1Mat(DVR this, Mat *M) {
  int ne; BPSGetNumEle(this->bps, &ne);
  int nq = this->nq;
  int n = (nq-2) * ne + ne -1;
  MatCreate(this->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode DVRSetSR1Mat(DVR this, Mat *M) {
  PetscErrorCode ierr;
  DVRInitR1Mat(this, M);
  
  for(int i = 0; i < this->num_basis; i++) 
    ierr = MatSetValue(*M, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRSetD2R1Mat(DVR this, Mat *M) {
  Mat LSR1;
  DVRSetD2R1LSMat(this, &LSR1);
  DVRLSMatToMat(this, LSR1, M);  
  MatDestroy(&LSR1);
  return 0;
}
PetscErrorCode DVRSetR2invR1Mat(DVR this, Mat *M) {
  PetscErrorCode ierr;
  DVRInitR1Mat(this, M);
  
  for(int i = 0; i < this->num_basis; i++) {
    PetscScalar x = this->xs_basis[i];
    ierr = MatSetValue(*M, i, i, 1.0/(x*x), INSERT_VALUES); CHKERRQ(ierr);
  }
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRSetENR1Mat(DVR this, int q, double a, Mat *M) {
  PetscErrorCode ierr;
  DVRInitR1Mat(this, M);
  
  for(int i = 0; i < this->num_basis; i++) {
    PetscScalar v;
    PartialCoulomb(q, a, this->xs_basis[i], &v);
    ierr = MatSetValue(*M, i, i, v, INSERT_VALUES); CHKERRQ(ierr);
  }
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}

// ------- R2Mat -----------
PetscErrorCode DVRSetEER2Mat(DVR this, int q, Mat *M) {
  PetscErrorCode ierr;
  Mat LS;
  ierr = DVRSetEER2LSMat(this, q, &LS); CHKERRQ(ierr);
  ierr = DVRR2LSMatToR2Mat(this, LS, M); CHKERRQ(ierr);
  return 0;
}

// -------- LSR1Mat --------
PetscErrorCode DVRInitR1LSMat(DVR this, Mat *M) {
  
  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);
  int n = nq * ne;
  MatCreate(this->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode DVRSetSR1LSMat(DVR this, Mat *M) {

  DVRInitR1LSMat(this, M);

  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);
  int n = nq * ne;

  for(int k = 0; k < n; k++)
    MatSetValue(*M, k, k, this->ws[k], INSERT_VALUES);

  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVRPrepareD2R1LSMat(DVR this) {

  Mat M;
  PetscErrorCode ierr;
  PetscScalar *ds_list;
  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);
  ds_list = (PetscScalar*)malloc(sizeof(PetscScalar)*nq*ne*nq);
  for(int i = 0; i < ne; i++)
    for(int m = 0; m < nq; m++)
      for(int mp = 0; mp < nq; mp++) {
	int idx = i*nq*nq + m*nq + mp;
	DerivLS(this->xs, this->ws, ne, nq, i, m, mp, &ds_list[idx]);
      }

  DVRInitR1LSMat(this, &M);

  for(int i = 0; i < ne; i++) 
    for(int ma = 0; ma < nq; ma++)
      for(int mb = 0; mb < nq; mb++) {
	int idx_a = i*nq + ma;
	int idx_b = i*nq + mb;
	PetscScalar v = 0.0;
	for(int k = 0; k < nq; k++) 
	  v -= ds_list[idx_a*nq+k] * ds_list[idx_b*nq+k] * this->ws[i*nq+k];
	MatSetValue(M, idx_a, idx_b, v, INSERT_VALUES);
      }

  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  free(ds_list);  
  this->D2_R1LSMat = M;
  return 0;  
}
PetscErrorCode DVRSetD2R1LSMat(DVR this, Mat *M){
  
  if(this->D2_R1LSMat == NULL) {
    PetscErrorCode ierr;
    ierr = DVRPrepareD2R1LSMat(this); CHKERRQ(ierr);
  }
  
  PetscErrorCode ierr;
  ierr = MatConvert(this->D2_R1LSMat, MATSAME, MAT_INITIAL_MATRIX, M); 
  CHKERRQ(ierr);  
  return 0;
  
}
PetscErrorCode DVRPrepareR2invR1LSMat(DVR this) {

  Mat M;
  PetscErrorCode ierr;
  DVRInitR1LSMat(this, &M);

  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);

  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) 
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w, v;
      x = this->xs[idx];
      w = this->ws[idx];
      if(i_ele == 0 && i_q == 0)
	v = 0.0;
      else
	v = 1.0/(x*x);

      ierr = MatSetValue(M, idx, idx, v*w, INSERT_VALUES); CHKERRQ(ierr);
      idx++;
    }    
  
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  this->R2_R1LSMat = M;

  return 0;  
}
PetscErrorCode DVRSetR2invR1LSMat(DVR this, Mat *M) {

  if(this->R2_R1LSMat == NULL) {
    PetscErrorCode ierr;
    ierr = DVRPrepareR2invR1LSMat(this); CHKERRQ(ierr);
  }
  
  PetscErrorCode ierr;
  ierr = MatConvert(this->R2_R1LSMat, MATSAME, MAT_INITIAL_MATRIX, M); 
  CHKERRQ(ierr);  
  return 0;
  
}
PetscErrorCode DVRSetENR1LSMat(DVR this, int q, double a, Mat *M) {

  PetscErrorCode ierr;

  DVRInitR1LSMat(this, M);

  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);

  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w, v;
      x = this->xs[idx];
      w = this->ws[idx];
      if(i_ele == 0 && i_q == 0)
	v = 0.0;
      else
	PartialCoulomb(q, a, x, &v);

      ierr = MatSetValue(*M, idx, idx, v*w, INSERT_VALUES); CHKERRQ(ierr);
      idx++;
    }    
  }  

  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  
  return 0;
}

PetscErrorCode DVRLSMatToMat(DVR this, Mat A, Mat *B) {

  MatMatMatMult(this->TT, A, this->T, MAT_INITIAL_MATRIX, PETSC_DEFAULT, B);
  return 0;  
  
}

// ---------- LSR2Mat ---------  
PetscErrorCode DVRInitR2LSMat(DVR this, Mat *M) {

  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);

  int n = pow(nq * ne, 2);
  MatCreate(this->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode DVRSetEER2LSMat(DVR this, int q, Mat *M) {

  PetscErrorCode ierr;

  ierr = DVRInitR2LSMat(this, M); CHKERRQ(ierr);

  Mat D2, R2;
  ierr = DVRSetD2R1LSMat(this, &D2);
  ierr = DVRSetR2invR1LSMat(this, &R2);
  
  int nq = this->nq;
  int ne; BPSGetNumEle(this->bps, &ne);
  PetscScalar rmax; BPSGetZMax(this->bps, &rmax);
  
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
    dgetrf_(&lnq, &lnq, vs, &lnq, ipiv, &info);
    dgetri_(&lnq, vs, &lnq, ipiv, work, &lnq, &info);

    for(int j_ele = 0; j_ele < ne; j_ele++) {
      for(int i = 0; i < nq; i++) 
	for(int j = 0; j < nq; j++) {
	  double ri = this->xs[i_ele*nq+i]; double wi = this->ws[i_ele*nq+i];
	  double rj = this->xs[j_ele*nq+j]; double wj = this->ws[j_ele*nq+j];
	  int idxA = (i_ele*nq+i)*nq*ne + (j_ele*nq+j);
	  double v = pow(ri, q) * pow(rj, q) / pow(rmax, 2*q+1);
	  if(i_ele == j_ele) 
	    if(i != 0 && j != 0)
	      v += (2*q+1) * vs[i*nq+j] / (ri*rj*sqrt(wi)*sqrt(wj));
	  ierr = MatSetValue(*M, idxA, idxA, v, INSERT_VALUES); CHKERRQ(ierr);	  
	}
    }
  }

  ierr = MatDestroy(&T); CHKERRQ(ierr);  
  PetscFree(vs); PetscFree(idx); PetscFree(work); PetscFree(ipiv);

  ierr = MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  
  ierr = MatAssemblyEnd(  *M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  
  
  return 0;
}

PetscErrorCode DVRR2LSMatToR2Mat(DVR this, Mat A, Mat *B) {

  PetscErrorCode ierr;

  if(this->T2 == NULL || this->T2T) 
    DVRPrepareT2(this);
  
  ierr = MatMatMatMult(this->T2T, A, this->T2, 
		       MAT_INITIAL_MATRIX, PETSC_DECIDE, B); CHKERRQ(ierr);  
  return 0;
}

