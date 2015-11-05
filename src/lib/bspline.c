#include <rescol/bspline.h>

// ---- External functions ----
int NumBSpline(int order, int num_ele) {
  return num_ele + 2*(order-1) - 2 - (order-1);
}
int HasNon0Value(int order, int i, int j) {
  return abs(i-j) < order;
}
PetscErrorCode CalcBSpline(int order, double* ts, int i, double x, double* y) {

  if(order < 1) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, 
	    "order must be positive integer\n");
  }

  if(x < ts[i] || ts[i+order] < x) {
    *y = 0.0; 
    return 0;
  }

  if(order == 1) {
    if(ts[i] <= x && x < ts[i+1])
      *y = 1.0;
    else
      *y = 0.0;
  } else {
    double ti    = ts[i];
    double ti1   = ts[i+1];
    double tidm1 = ts[i+order-1];
    double tid   = ts[i+order];
    double bs0; CalcBSpline(order-1, ts, i, x, &bs0);
    double bs1; CalcBSpline(order-1, ts, i+1, x, &bs1);
    
    double acc = 0.0;
    if(fabs(bs0) > 0.000000001)
      acc += (x - ti) / (tidm1 - ti) * bs0;
    if(fabs(bs1) > 0.000000001)
      acc += (tid - x) / (tid - ti1) * bs1;
    *y = acc;
  }
  return 0;
}
PetscErrorCode CalcDerivBSpline(int order, double* ts, int i, double x, double* y) {
  
  if(order < 1) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, 
	    "order must be positive integer\n");
  }
  
  if(order == 1) 
    *y = 0.0;
  else {

    int k = order;
    double bs0;
    CalcBSpline(k-1, ts, i, x, &bs0);
    double bs1;
    CalcBSpline(k-1, ts, i+1, x, &bs1);
    double eps = 0.000000001;
    double acc = 0.0;
    if(fabs(bs0) > eps)
      acc += (k-1)/(ts[i+k-1]-ts[i]) * bs0;
    if(fabs(bs1) > eps)
      acc -= (k-1)/(ts[i+k]-ts[i+1]) * bs1;
    *y = acc;
  }
  return 0;
}
PetscErrorCode Non0QuadIndex(int a, int c, int k, int nq, int* i0, int* i1) {

  int a0_idx, a1_idx;
  a0_idx = a < c ? a + 2 : c + 2;
  a1_idx = a < c ? c + 2 : a + 2;

  *i0 = a1_idx*k - k*k;
  *i1 = a0_idx*k;
  if(*i0<0)
    *i0 = 0;
  if(*i1>nq)
    *i1 = nq;

  /*
  *i0 = a<c ? (c-k+2)*k : (a-k+2)*k;
  if(*i0<0)
    *i0 = 0;
  *i1 = a<c ? (a+k-1)*k : (c+k-1)*k;
  if(*i1>nq)
    *i1=nq;
    */
  return 0;
}

// ---- Basic Methods ----
PetscErrorCode BSSCreate(BSS *bss, int order, BPS bps, Scaler scaler, 
			 MPI_Comm comm) {
  PetscErrorCode ierr;
  int i, ib, ie, iq;
  BSS _bss;

  PetscNew(&_bss);
  *bss = NULL;

  PetscScalar *zs; PetscInt num_zs;
  BPSGetZs(bps, &zs, &num_zs);
  
  // data num
  _bss->comm = comm;
  _bss->order = order;
  _bss->bps = bps;
  if(scaler == NULL) {
    Scaler scaler_none;
    ScalerCreateNone(&scaler_none, comm);
    _bss->scaler = scaler_none;
  } else
    _bss->scaler = scaler;

  BPSGetNumEle(bps, &_bss->num_ele);
  _bss->num_basis = NumBSpline(order, num_zs-1);
  //  BPSGetZMax(bps, &_bss->rmax);

  // copy ts and zs
  PetscMalloc1(num_zs+2*order-2, &_bss->ts);

  for(i = 0; i < order-1; i++) {
    _bss->ts[i] = zs[0];
    _bss->ts[order-1+num_zs+i] = zs[num_zs-1];
  }
  for(i = 0; i < num_zs; i++) 
    _bss->ts[i+order-1] = zs[i];

  // calculate appreciate quadrature points
  int n_xs = _bss->num_ele * _bss->order;
  PetscMalloc1(_bss->num_basis, &_bss->b_idx_list);
  PetscMalloc1(n_xs, &_bss->xs);
  PetscMalloc1(n_xs, &_bss->ws);
  PetscMalloc1(n_xs, &_bss->qrs);
  PetscMalloc1(n_xs, &_bss->Rrs);
  int num = sizeof(PetscScalar)*(n_xs)*(_bss->num_basis);
  PetscMalloc1(num, &_bss->vals);
  PetscMalloc1(num, &_bss->derivs);

  for(ib = 0; ib < _bss->num_basis; ib++)
    _bss->b_idx_list[ib] = ib + 1;
  for(ie = 0; ie < _bss->num_ele; ie++) {
    PetscScalar a, b; a = zs[ie]; b = zs[ie+1];
    for(iq = 0; iq < order; iq++) {
      PetscScalar x, w;
      int ix = ie*order+iq;
      LegGauss(order, iq, &x, &w);
      x = (b+a)/2.0  + (b-a)/2.0 * x; w = (b-a)/2.0*w;
      _bss->xs[ix] = x;
      _bss->ws[ix] = w;
      
      for(ib = 0; ib < _bss->num_basis; ib++) {
	PetscReal y;
	PetscReal dy;
	CalcBSpline(     order, _bss->ts, _bss->b_idx_list[ib], x, &y);
	CalcDerivBSpline(order, _bss->ts, _bss->b_idx_list[ib], x, &dy);
	int iy = ib*(_bss->num_ele*order) + ie*order + iq;
	_bss->vals[iy] = y;  
	_bss->derivs[iy] = dy;
      }
    }
  }
  
  ierr = ScalerSetQr(_bss->scaler, _bss->xs, n_xs, _bss->qrs); CHKERRQ(ierr);
  ierr = ScalerSetRr(_bss->scaler, _bss->xs, n_xs, _bss->Rrs); CHKERRQ(ierr);

  PetscFree(zs);
  *bss = _bss;
  return 0;
}
PetscErrorCode BSSCreateFromOptions(BSS *bss, MPI_Comm comm) {
  PetscBool find;
  PetscInt order;
  BPS bps;
  Scaler scaler;
  PetscErrorCode ierr;

  order = 2;
  ierr = PetscOptionsGetInt(NULL, "-bss_order", &order, &find); CHKERRQ(ierr);
  ierr = BPSCreate(&bps, comm); BPSSetFromOptions(bps); CHKERRQ(ierr);
  ierr = ScalerCreateFromOptions(&scaler, comm); CHKERRQ(ierr);
  ierr = BSSCreate(bss, order, bps, scaler, comm);  CHKERRQ(ierr);

  return 0;
 }
PetscErrorCode BSSDestroy(BSS *bss) {
  PetscErrorCode ierr;
  BSS this = *bss;
  ierr = BPSDestroy(&this->bps); CHKERRQ(ierr);
  ierr = ScalerDestroy(&this->scaler); CHKERRQ(ierr);
  ierr = PetscFree(this->b_idx_list); CHKERRQ(ierr);
  ierr = PetscFree(this->ts);  CHKERRQ(ierr);
  ierr = PetscFree(this->xs); CHKERRQ(ierr);
  ierr = PetscFree(this->ws); CHKERRQ(ierr);
  ierr = PetscFree(this->vals); CHKERRQ(ierr);
  ierr = PetscFree(this->derivs); CHKERRQ(ierr);
  ierr = PetscFree(*bss); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode BSSFPrintf(BSS this, FILE* file, int lvl) {

  MPI_Comm comm = this->comm;

  if(lvl != 0) {
    SETERRQ(comm, PETSC_ERR_ARG_OUTOFRANGE, 
	    "now only lvl=0 is supported.");
  }

  PetscFPrintf(comm, file, "===== Begin B-Spline =====\n");
  PetscFPrintf(comm, file, "order: %d\n", this->order);
  PetscFPrintf(comm, file, "num_basis: %d\n", this->num_basis);
  BPSFPrintf(this->bps, file, lvl);
  PetscFPrintf(comm, file, "===== End B-Spline =====\n");
  return 0;
}
PetscErrorCode BSSBasisPsi(BSS this, int i, PetscReal x, PetscReal *y) {

  PetscReal z;
  CalcBSpline(this->order, 
	      this->ts, 
	      this->b_idx_list[i], 
	      x, 
	      &z);
  *y = z;
  return 0;
}
PetscErrorCode BSSDerivBasisPsi(BSS this, int i, PetscReal x, PetscReal *y) {
  PetscReal z;
  CalcDerivBSpline(this->order, 
		   this->ts, 
		   this->b_idx_list[i], 
		   x, 
		   &z);
  *y = z;
  return 0;
}

// ---- Accessor ----
PetscErrorCode BSSGetSize(BSS this, int *n) {
 
  *n = this->num_basis;
  return 0;
}

// ---- Matrix -----
PetscErrorCode BSSInitR1Mat(BSS this, Mat *M) {
  int nb = this->num_basis;
  MatCreate(this->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nb, nb);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode BSSInitR2Mat(BSS this, Mat *M) {
  int nb = this->num_basis;
  MatCreate(this->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nb*nb, nb*nb);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode BSSCalcSR1Mat(BSS this, Mat S, InsertMode mode) {
  int i, j, k;
  int nb = this->num_basis;
  int ne = this->num_ele;
  int nq = this->order;
  PetscErrorCode ierr;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(this->order, this->b_idx_list[i], this->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++) 
	  v += this->vals[k+i*(ne*nq)] * this->vals[k+j*(ne*nq)] * this->ws[k] * this->qrs[k];
	ierr = MatSetValue(S, i, j, v, mode); CHKERRQ(ierr);
      }
    }
  return 0;
}
PetscErrorCode BSSCalcR2invR1Mat(BSS this, Mat M, InsertMode mode) {

  int i, j, k;
  int nb = this->num_basis;
  int ne = this->num_ele;
  int nq = this->order;
  PetscErrorCode ierr;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(this->order, this->b_idx_list[i], this->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++) {
	  PetscScalar x = this->Rrs[k];
	  PetscScalar w = this->ws[k];
	  PetscScalar q = this->qrs[k];
	  PetscScalar fi = this->vals[k+i*ne*nq];
	  PetscScalar fj = this->vals[k+j*ne*nq];
	  v += fi*fj*w*q/(x*x);
	}
	ierr = MatSetValue(M, i, j, v, mode); CHKERRQ(ierr);
      }
    }
  return 0;
}
PetscErrorCode BSSCalcD2R1Mat(BSS this, Mat D, InsertMode mode) {
  int i, j, k;
  int nb = this->num_basis;
  int ne = this->num_ele;
  int nq = this->order;
  PetscErrorCode ierr;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(this->order, this->b_idx_list[i], this->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++) 
	  v += this->derivs[k+i*(ne*nq)] * this->derivs[k+j*(ne*nq)] * this->ws[k] / this->qrs[k];
	ierr = MatSetValue(D, i, j, -v, mode); CHKERRQ(ierr);
      }
    }
  return 0;
}
PetscErrorCode BSSCalcENR1Mat(BSS this, int q, PetscReal a, Mat V, InsertMode mode) {
  int order = this->order;
  int nb = this->num_basis;
  int ne = this->num_ele;
  int nq = order*ne;
  PetscErrorCode ierr;

  PetscScalar *vs; 
  PetscMalloc1(nq, &vs);
  for(int k = 0; k < nq; k++) {
    PetscScalar g = this->xs[k] > a ? this->Rrs[k] : a;
    PetscScalar s = this->xs[k] < a ? this->Rrs[k] : a;
    vs[k] = pow(s/g, q)/g * this->ws[k] * this->qrs[k];
  }

  for(int i = 0; i < nb; i++) {
    int j0 = i-order+1; j0 = j0 < 0 ? 0: j0;
    for(int j = j0; j <= i; j++) {
      int k0, k1;  
      Non0QuadIndex(j, i, order, nq, &k0, &k1);
      PetscScalar v = 0.0;
      for(int k = k0; k < k1; k++) 
	v += this->vals[k+i*nq] * this->vals[k+j*nq] * vs[k];
      ierr = MatSetValue(V, i, j, v, mode); CHKERRQ(ierr);
      if(i!=j)
	ierr = MatSetValue(V, j, i, v, mode); CHKERRQ(ierr);
    }
  }

  PetscFree(vs);
  return 0;
}
PetscErrorCode BSSCalcEER2Mat(BSS this, int q, Mat V, InsertMode mode) {

  int k = this->order;
  int nb = this->num_basis;
  int nq = this->order * this->num_ele;
  PetscErrorCode ierr;

  PetscScalar *sg_ij; 
  ierr = PetscMalloc1(nq*nq, &sg_ij); CHKERRQ(ierr);
  for(int i = 0; i < nq; i++)
    for(int j = 0; j < nq; j++) {
      PetscScalar g = this->xs[i]>this->xs[j] ? this->Rrs[i] : this->Rrs[j];
      PetscScalar s = this->xs[i]>this->xs[j] ? this->Rrs[j] : this->Rrs[i];
      sg_ij[j+nq*i] = pow(s/g, q)/g;      
    }

  PetscReal *wvv; ierr = PetscMalloc1(nb*nb*nq, &wvv); CHKERRQ(ierr);
  int *i0s; ierr = PetscMalloc1(nb*nb, &i0s); CHKERRQ(ierr);
  int *i1s; ierr = PetscMalloc1(nb*nb, &i1s); CHKERRQ(ierr);
  for(int a = 0; a < nb; a++){
    int c0 = a-k+1; c0 = c0<0?0:c0;
    int c1 = a+k  ; c1 = c1>nb?nb:c1;
    for(int c = c0; c < c1; c++) {
      int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
      i0s[a*nb+c] = i0; i1s[a*nb+c] = i1;
      for(int n = 0; n < nq; n++) 
	wvv[a*nb*nq+c*nq+n]= this->ws[n] * this->vals[a*nq+n] * this->vals[c*nq+n];
    }
  }

  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1; c0 = c0<0?0:c0;
    for(int c = c0; c <= a; c++) {
      for(int b = 0; b < nb; b++) {
	int d0 = b-k+1; d0 = d0<0?0:d0;
	for(int d = d0; d <= b; d++) {
	  double v = 0.0;
	  for(int i = i0s[a*nb+c]; i < i1s[a*nb+c]; i++) {
	    double aci = wvv[a*nb*nq + c*nq + i];
	    for(int j = i0s[b*nb+d]; j < i1s[b*nb+d]; j++)
	      v += aci * sg_ij[i*nq+j] * wvv[b*nb*nq + d*nq + j];
	  }
	  
	  if(a > c && b > d) {
	    ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(V, a*nb+d, c*nb+b, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(V, c*nb+b, a*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(V, c*nb+d, a*nb+b, v, mode); CHKERRQ(ierr);
	  } else if(a > c && b == d) {
	    ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(V, c*nb+b, a*nb+d, v, mode); CHKERRQ(ierr);
	  } else if(a == c && b > d) {
	    ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(V, a*nb+d, c*nb+b, v, mode); CHKERRQ(ierr);   
	  } else if(a == c && b == d) {
	    ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	  }
	}
      }
    }
  }

  PetscFree(sg_ij); PetscFree(wvv); PetscFree(i0s); PetscFree(i1s);
  return 0;
}
PetscErrorCode BSSCalcEER2Mat_ver1(BSS this, int q, Mat V, InsertMode mode) {

  int nb = this->num_basis;
  int nq = this->order * this->num_ele;
  PetscErrorCode ierr;

  double *sg_ij;
  sg_ij = (double*)malloc(sizeof(double)*nq*nq);
  for(int i = 0; i < nq; i++)
    for(int j = 0; j < nq; j++) {
      double v; PartialCoulomb(q, this->xs[i], this->xs[j], &v);
      sg_ij[j+nq*i] = v;
    }

  for(int a = 0; a < nb; a++) {
    // int c0 = a-k+1; c0 = c0<0?0:c0;
    // int c1 = a+k  ; c1 = c1>nb?nb:c1;
    int c0 = 0; int c1 = nb;
    for(int c = c0; c < c1; c++) {
      for(int b = 0; b < nb; b++) {
	//int d0 = b-k+1; d0 = d0<0?0:d0;
	//int d1 = b+k;   d1 = d1>nb?nb:d1;
	int d0 = 0; int d1 = nb;
	for(int d = d0; d < d1; d++) {
	  double v = 0.0;
	  for(int n = 0; n < nq; n++) {
	    for(int m = 0; m < nq; m++){
	      v += this->ws[n] * this->vals[a*nq+m] * this->vals[c*nq+m] 
		* this->ws[m] * this->vals[b*nq+n] * this->vals[d*nq+n] 
		* sg_ij[n*nq+m];
	    }
	  }
	  ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	}
      }
    }
  }
  return 0;
}
PetscErrorCode BSSSetSR1Mat(BSS this, Mat *S) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, S); CHKERRQ(ierr);
  ierr = BSSCalcSR1Mat(this, *S, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*S, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSSetR2invR1Mat(BSS this, Mat *M) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, M); CHKERRQ(ierr);
  ierr = BSSCalcR2invR1Mat(this, *M, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;  
}
PetscErrorCode BSSSetD2R1Mat(BSS this, Mat* D) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, D); CHKERRQ(ierr);
  ierr = BSSCalcD2R1Mat(this, *D, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*D, MAT_FINAL_ASSEMBLY);
  return 0;  
}
PetscErrorCode BSSSetENR1Mat(BSS this, int q, PetscScalar a, Mat *D) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, D); CHKERRQ(ierr);
  ierr = BSSCalcENR1Mat(this, q, a, *D, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*D, MAT_FINAL_ASSEMBLY);
  return 0;  
}
PetscErrorCode BSSSetPotR1Mat(BSS this, POT pot, Mat *M) {
  PetscErrorCode ierr;
  int nb = this->num_basis;
  int ne = this->num_ele;
  int nq = this->order;
  BSSInitR1Mat(this, M);

  PetscScalar *vs;
  PetscMalloc1(ne*nq, &vs);

  for(int k = 0; k < ne*nq; k++) 
    vs[k] = POTCalc(pot, this->xs[k]);

  for(int i = 0; i < nb; i++)
    for(int j = 0; j < nb; j++) {
      if(HasNon0Value(this->order, this->b_idx_list[i], this->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(int k = 0; k < ne*nq; k++) {
	  PetscReal w = this->ws[k];
	  PetscScalar v1 = vs[k];
	  v += this->vals[k+i*(ne*nq)] * this->vals[k+j*(ne*nq)] 
	    * w * v1;
	}
	ierr = MatSetValue(*M, i, j, v, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  PetscFree(vs);
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSSetEER2Mat(BSS this, int q, Mat *V) {
  PetscErrorCode ierr;
  ierr = BSSInitR2Mat(this, V); CHKERRQ(ierr);
  ierr = BSSCalcEER2Mat(this, q, *V, INSERT_VALUES ); CHKERRQ(ierr);
  MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);
  return 0;
}

