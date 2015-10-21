#include "bspline.h"


// External functions
int NumBSpline(int order, int num_ele) {
  return num_ele + 2*(order-1) - 2 - (order-1);
}

int HasNon0Value(int order, int i, int j) {
  return abs(i-j) < order;
}

PetscErrorCode PartialCoulomb(int q, double r1, double r2, double *y) {

  if(q < 0) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, 
	    "q must be non negative integer");
  }

  double g = r1>r2 ? r1 : r2;
  double s = r1>r2 ? r2 : r1;
  *y = pow(s/g, q)/g;
  return 0;
}

PetscScalar LegGauss(int n, int i, PetscScalar* x, PetscScalar* w) {
  PetscScalar xs[45] = {
    0.0, -0.5773502691896257, 0.5773502691896257, -0.7745966692414834, 0.0, 0.7745966692414834, -0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526, -0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664, -0.932469514203152, -0.6612093864662645, -0.23861918608319693, 0.23861918608319693, 0.6612093864662645, 0.932469514203152, -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585, -0.9602898564975362, -0.7966664774136267, -0.525532409916329, -0.18343464249564984, 0.18343464249564984, 0.525532409916329, 0.7966664774136267, 0.9602898564975362, -0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
  PetscScalar ws[45] = {
    2.0, 1.0, 1.0, 0.5555555555555555, 0.888888888888889, 0.5555555555555555, 0.34785484513745396, 0.6521451548625462, 0.6521451548625462, 0.34785484513745396, 0.236926885056189, 0.4786286704993665, 0.5688888888888888, 0.4786286704993665, 0.236926885056189, 0.17132449237916944, 0.36076157304813883, 0.4679139345726918, 0.4679139345726918, 0.36076157304813883, 0.17132449237916944, 0.12948496616886826, 0.2797053914892774, 0.38183005050511937, 0.41795918367347007, 0.38183005050511937, 0.2797053914892774, 0.12948496616886826, 0.10122853629037527, 0.22238103445337512, 0.3137066458778874, 0.362683783378362, 0.362683783378362, 0.3137066458778874, 0.22238103445337512, 0.10122853629037527, 0.08127438836157427, 0.18064816069485748, 0.2606106964029357, 0.31234707704000275, 0.3302393550012596, 0.31234707704000275, 0.2606106964029357, 0.18064816069485748, 0.0812743883615742};

  *x = xs[n * (n-1)/2 + i];
  *w = ws[n * (n-1)/2 + i];
  return 0;
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

PetscErrorCode CreateLinKnots(int num, double zmax, double *zs[]) {
  double dz = zmax / (num-1);
  *zs = (double*)malloc(sizeof(double)*num);
  for(int i = 0; i < num; i++) 
    (*zs)[i] = i * dz;
  return 0;
}

PetscErrorCode CreateExpKnots(int num, double zmax, double gamma, double *zs[]){
  *zs = (double*)malloc(sizeof(double)*num);
  for(int n = 0; n < num; n++) {
    (*zs)[n] = zmax * (exp(gamma*n/(num-1)) - 1.0) / (exp(gamma) - 1.0);
  }
  return 0;
}

PetscErrorCode Non0QuadIndex(int a, int c, int k, int nq, int* i0, int* i1) {
  *i0 = a<c ? (c-k+2)*k : (a-k+2)*k;
  if(*i0<0)
    *i0 = 0;
  *i1 = a<c ? (a+k-1)*k : (c+k-1)*k;
  if(*i1>nq)
    *i1=nq;

  return 0;
}

// Methods
PetscErrorCode BSSCreate(BSS *bss, int order, double*zs, int num_zs) {
  int i, ib, ie, iq;
  BSS _bss;

  _bss = (BSS)malloc(sizeof(struct _p_BSS));
  *bss = NULL;
  
  // data num
  _bss->order = order;
  _bss->num_ele = num_zs-1;
  _bss->num_basis = NumBSpline(order, num_zs-1);
  _bss->rmax = zs[num_zs-1];
  strcpy(_bss->knots_type, "unknown");

  // copy ts and zs
  _bss->zs = (PetscScalar*)malloc(sizeof(PetscScalar)*(num_zs));
  _bss->ts = (PetscScalar*)malloc(sizeof(PetscScalar)*(num_zs+2*order-2));

  for(i = 0; i < order-1; i++) {
    _bss->ts[i] = zs[0];
    _bss->ts[order-1+num_zs+i] = zs[num_zs-1];
  }
  for(i = 0; i < num_zs; i++) {
    _bss->zs[i] = zs[i];
    _bss->ts[i+order-1] = zs[i];
  }

  // calculate appreciate quadrature points
  int n_xs = _bss->num_ele * _bss->order;
  _bss->b_idx_list = (int*)malloc(sizeof(int)*(_bss->num_basis));
  _bss->xs = (PetscScalar*)malloc(sizeof(PetscScalar)*n_xs);
  _bss->ws = (PetscScalar*)malloc(sizeof(PetscScalar)*n_xs);
  int num = sizeof(PetscScalar)*(n_xs)*(_bss->num_basis);
  _bss->vals = (PetscScalar*)malloc(num);
  _bss->derivs = (PetscScalar*)malloc(num);

  for(ib = 0; ib < _bss->num_basis; ib++)
    _bss->b_idx_list[ib] = ib + 1;
  for(ie = 0; ie < _bss->num_ele; ie++) {
    PetscScalar a, b; a = _bss->zs[ie]; b = _bss->zs[ie+1];
    for(iq = 0; iq < order; iq++) {
      PetscScalar x, w;
      int ix = ie*order+iq;
      LegGauss(order, iq, &x, &w);
      x = (b+a)/2.0  + (b-a)/2.0 * x; w = (b-a)/2.0*w;
      _bss->xs[ix] = x;
      _bss->ws[ix] = w;
      
      for(ib = 0; ib < _bss->num_basis; ib++) {
	PetscScalar y;
	PetscScalar dy;
	CalcBSpline(     order, _bss->ts, _bss->b_idx_list[ib], x, &y);
	CalcDerivBSpline(order, _bss->ts, _bss->b_idx_list[ib], x, &dy);
	int iy = ib*(_bss->num_ele*order) + ie*order + iq;
	_bss->vals[iy] = y; _bss->derivs[iy] = dy;
      }
    }
  }

  *bss = _bss;
  return 0;
}

PetscErrorCode BSSCreateFromOptions(BSS *bss, MPI_Comm comm) {
  PetscBool find;
  PetscReal rmax;
  PetscInt order, num;
  PetscErrorCode ierr;
  char knots[10] = "line";

  order = 2;
  ierr = PetscOptionsGetInt(NULL, "-bss_order", &order, &find); CHKERRQ(ierr);
  rmax = 20.0;
  ierr = PetscOptionsGetReal(NULL, "-bss_rmax", &rmax, &find); CHKERRQ(ierr);
  num = 21;
  ierr = PetscOptionsGetInt(NULL, "-bss_knots_num", &num, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-bss_knots_type", knots, 10, &find); 
  CHKERRQ(ierr);
  
  double *zs; 
  char knots_type[10];
  if(strcmp(knots, "line") == 0) {
    ierr = CreateLinKnots(num, rmax, &zs); CHKERRQ(ierr);
    strcpy(knots_type, "line");
  } else if(strcmp(knots, "exp") == 0) {
    ierr = CreateExpKnots(num, rmax, 5.0, &zs); CHKERRQ(ierr);
    strcpy(knots_type, "exp");
  } else {
    SETERRQ(comm, 1, "bss_knots_type must be line or exp."); }

  ierr = BSSCreate(bss, order, zs, num);  CHKERRQ(ierr);
  strcpy((*bss)->knots_type, knots_type);
  return 0;
 }

PetscErrorCode BSSDestroy(BSS *bss) {
   BSS this = *bss;
   free(this->ts); 
   free(this->xs);
   free(this->vals);
   free(this->derivs);
   return 0;
}

PetscErrorCode BSSFPrintf(BSS this, MPI_Comm comm, FILE* file, int lvl) {

  if(lvl != 0) {
    SETERRQ(comm, PETSC_ERR_ARG_OUTOFRANGE, 
	    "now only lvl=0 is supported.");
  }

  PetscFPrintf(comm, file, "SUMMARY: B-Spline set\n");
  PetscFPrintf(comm, file, "order: %d\n", this->order);
  PetscFPrintf(comm, file, "num_ele: %d\n", this->num_ele);
  PetscFPrintf(comm, file, "num_basis: %d\n", this->num_basis);
  PetscFPrintf(comm, file, "knots_type: %s\n", this->knots_type);
  PetscFPrintf(comm, file, "rmax: %f\n", this->rmax);
  return 0;
}

PetscErrorCode BSSBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y) {
  PetscScalar z;
  CalcBSpline(this->order, 
	      this->ts, 
	      this->b_idx_list[i], 
	      x, 
	      &z);
  *y = z;
  return 0;
}

PetscErrorCode BSSDerivBasisPsi(BSS this, int i, PetscScalar x, PetscScalar *y) {
  PetscScalar z;
  CalcDerivBSpline(this->order, 
		   this->ts, 
		   this->b_idx_list[i], 
		   x, 
		   &z);
  *y = z;
  return 0;
}

PetscErrorCode BSSInitR1Mat(BSS this, MPI_Comm comm, Mat *M) {
  int nb = this->num_basis;
  MatCreate(comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nb, nb);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;
}

PetscErrorCode BSSInitR2Mat(BSS this, MPI_Comm comm, Mat *M) {
  int nb = this->num_basis;
  MatCreate(comm, M);
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
	  v += this->vals[k+i*(ne*nq)] * this->vals[k+j*(ne*nq)] * this->ws[k];
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
	  PetscScalar x = this->xs[k];
	  PetscScalar w = this->ws[k];
	  PetscScalar fi = this->vals[k+i*ne*nq];
	  PetscScalar fj = this->vals[k+j*ne*nq];
	  v += fi*fj*w/(x*x);
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
	  v += this->derivs[k+i*(ne*nq)] * this->derivs[k+j*(ne*nq)] * this->ws[k];
	ierr = MatSetValue(D, i, j, -v, mode); CHKERRQ(ierr);
      }
    }
  return 0;
}

PetscErrorCode BSSCalcENR1Mat(BSS this, int q, double a, Mat V, InsertMode mode) {
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
	  double r = this->xs[k];
	  double w = this->ws[k];
	  double v1;
	  PartialCoulomb(q, a, r, &v1);
	  v += this->vals[k+i*(ne*nq)] * this->vals[k+j*(ne*nq)] 
	    * w * v1;
	}
	ierr = MatSetValue(V, i, j, v, mode); CHKERRQ(ierr);
      }
    }
  return 0;
}
PetscErrorCode BSSCalcEER2Mat(BSS this, int q, Mat V, InsertMode mode) {

  int k = this->order;
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

  double *wvv = (double*)malloc(sizeof(double)*nb*nb*nq);
  int *i0s = (int*)malloc(sizeof(int)*nb*nb);
  int *i1s = (int*)malloc(sizeof(int)*nb*nb);
  for(int a = 0; a < nb; a++){
    int c0 = a-k+1; c0 = c0<0?0:c0;
    int c1 = a+k  ; c1 = c1>nb?nb:c1;
    for(int c = c0; c < c1; c++) {
      //int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
      //i0s[a*nb+c] = i0; i1s[a*nb+c] = i1;
      i0s[a*nb+c] = 0; i1s[a*nb+c] = nq;
      for(int n = 0; n < nq; n++) 
	wvv[a*nb*nq+c*nq+n]= this->ws[n] * this->vals[a*nq+n] * this->vals[c*nq+n];
    }
  }

  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1; c0 = c0<0?0:c0;
    int c1 = a+k  ; c1 = c1>nb?nb:c1;
    for(int c = c0; c < c1; c++) {
      for(int b = 0; b < nb; b++) {
	int d0 = b-k+1; d0 = d0<0?0:d0;
	int d1 = b+k;   d1 = d1>nb?nb:d1;
	for(int d = d0; d < d1; d++) {
	  double v = 0.0;
	  for(int n = i0s[a*nb+c]; n < i1s[a*nb+c]; n++) 
	    for(int m = i0s[b*nb+d]; m < i1s[b*nb+d]; m++)
	      v += wvv[a*nb*nq+c*nq+m] * wvv[b*nb*nq+d*nq+n]* sg_ij[n*nq+m];
	  ierr = MatSetValue(V, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	}
      }
    }
  }

  free(sg_ij);
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

PetscErrorCode BSSSetSR1Mat(BSS this, MPI_Comm comm, Mat *S) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, comm, S); CHKERRQ(ierr);
  ierr = BSSCalcSR1Mat(this, *S, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*S, MAT_FINAL_ASSEMBLY);
  return 0;
}

PetscErrorCode BSSSetR2invR1Mat(BSS this, MPI_Comm comm, Mat *M) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, comm, M); CHKERRQ(ierr);
  ierr = BSSCalcR2invR1Mat(this, *M, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;  
}

PetscErrorCode BSSSetD2R1Mat(BSS this, MPI_Comm comm, Mat* D) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, comm, D); CHKERRQ(ierr);
  ierr = BSSCalcD2R1Mat(this, *D, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*D, MAT_FINAL_ASSEMBLY);
  return 0;  
}

PetscErrorCode BSSSetENR1Mat(BSS this, int q, double a, MPI_Comm comm, Mat *D) {
  PetscErrorCode ierr;
  ierr = BSSInitR1Mat(this, comm, D); CHKERRQ(ierr);
  ierr = BSSCalcENR1Mat(this, q, a, *D, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*D, MAT_FINAL_ASSEMBLY);
  return 0;  
}

PetscErrorCode BSSSetEER2Mat(BSS this, int q, MPI_Comm comm, Mat *V) {
  PetscErrorCode ierr;
  ierr = BSSInitR2Mat(this, comm, V); CHKERRQ(ierr);
  ierr = BSSCalcEER2Mat(this, q, *V, INSERT_VALUES ); CHKERRQ(ierr);
  MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);
  return 0;
}

PetscErrorCode BSSSetUR1R2Mat(BSS this, Mat *U) {
  /*
    {<B_m | 1/r | rho_A>}_mA
    rho_A(r) = B_i(r)B_j(r)
   */
  
  return 0;
  

}

PetscErrorCode BSSSetEER2MatGreen(BSS this, int q, MPI_Comm, Mat *V) {
  return 0;
}
