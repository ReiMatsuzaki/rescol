#include "../include/bspline.h"

// ---- External functions ----
int NumBSpline(int order, int num_ele) {
  return num_ele + 2*(order-1) - 2 - (order-1);
}
int HasNon0Value(int order, int i, int j) {
  return abs(i-j) < order;
}
PetscErrorCode CalcBSpline(int k, double* ts_r, PetscScalar* ts, 
			   int i, double x_r,   PetscScalar x, PetscScalar* y,
			   PetscBool *zeroq) {

  if(k < 1) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, 
	    "order must be positive integer\n");
  }

  PetscBool _zeroq;

  if(x_r < ts_r[i] || ts_r[i+k] < x_r) {
    *y = 0.0; 
    if(zeroq != NULL)
      *zeroq = PETSC_TRUE;
    return 0;
  }

  if(k == 1) {
    if(ts_r[i] <= x_r && x_r < ts_r[i+1]) {
      *y = 1.0;
      _zeroq = PETSC_FALSE;
    } else {
      *y = 0.0;
      _zeroq = PETSC_TRUE;
    }
  } else {
    PetscScalar ti    = ts[i];
    PetscScalar ti1   = ts[i+1];
    PetscScalar tidm1 = ts[i+k-1];
    PetscScalar tid   = ts[i+k];
    PetscBool zeroq0, zeroq1;
    PetscScalar bs0; CalcBSpline(k-1, ts_r, ts, i,   x_r, x, &bs0, &zeroq0);
    PetscScalar bs1; CalcBSpline(k-1, ts_r, ts, i+1, x_r, x, &bs1, &zeroq1);
    
    _zeroq = PETSC_TRUE;
    PetscScalar acc = 0.0;
    //if(ScalarAbs(bs0) > 0.000000001 && !zeroq0) {
    if(!zeroq0) {
      PetscScalar x_minus_ti = x - ti;
      if(ScalarAbs(x_minus_ti) > 0.00000000001) {
	acc += x_minus_ti / (tidm1 - ti) * bs0;
      } else {
	zeroq0 = PETSC_TRUE;
      }
    }
    //    if(ScalarAbs(bs1) > 0.000000001 && !zeroq1) {
    if( !zeroq1) {
      if(ScalarAbs(tid-x) > 0.00000000001) {
	acc += (tid - x) / (tid - ti1) * bs1;
      } else {
	zeroq1 = PETSC_TRUE;
      }
    }
    *y = acc;
    _zeroq = zeroq1 && zeroq0;
  }
  if(zeroq != NULL)
    *zeroq = _zeroq;
  return 0;
}
PetscErrorCode CalcDerivBSpline(int order, PetscReal* ts_r, PetscScalar *ts, 
				int i, double x_r, PetscScalar x, PetscScalar* y,
				PetscBool *zeroq) {
  
  if(order < 1) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, 
	    "order must be positive integer\n");
  }
  
  PetscBool _zeroq;  

  if(order == 1) {
    *y = 0.0;
    _zeroq = PETSC_TRUE;
  } else {

    int k = order;
    PetscScalar bs0;
    PetscBool zeroq0;
    CalcBSpline(k-1, ts_r, ts, i, x_r, x, &bs0, &zeroq0);
    PetscScalar bs1;
    PetscBool zeroq1;
    CalcBSpline(k-1, ts_r, ts, i+1, x_r, x, &bs1, &zeroq1);
    double eps = 0.000000001;
    PetscScalar acc = 0.0;
    if(ScalarAbs(bs0) > eps && !zeroq0)
    acc += (k-1)/(ts[i+k-1]-ts[i]) * bs0;
    if(ScalarAbs(bs1) > eps && !zeroq1)
      acc -= (k-1)/(ts[i+k]-ts[i+1]) * bs1;
    *y = acc;
    _zeroq = (zeroq0 && zeroq1);
  }
  if(zeroq != NULL) 
    *zeroq = _zeroq;
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
PetscErrorCode BSSCreate(MPI_Comm comm, BSS *p_self) {
  
  BSS self;
  PetscNew(&self);
  self->comm = comm;
  self->order = 0;
  self->bps = NULL;
  self->c_scaling = NULL;

  self->num_ele = 0;
  self->num_basis  = 0;
  self->b_idx_list = NULL;
  self->ts_r = NULL;
  self->ts_s = NULL;
  self->xs = NULL;
  self->ws = NULL;
  self->qrs = NULL;
  self->Rrs = NULL;

  self->vals = NULL;
  self->derivs = NULL;

  self->set_knots = PETSC_FALSE;
  self->setup = PETSC_FALSE;
    
  *p_self = self;
  return 0;
}
PetscErrorCode BSSDestroy(BSS *p_self) {
  PetscErrorCode ierr;
  ierr = BSSCheck(*p_self); CHKERRQ(ierr);

  BSS self = *p_self;
  ierr = BPSDestroy(&self->bps); CHKERRQ(ierr);
  ierr = CScalingDestroy(&self->c_scaling); CHKERRQ(ierr);
  ierr = PetscFree(self->b_idx_list); CHKERRQ(ierr);
  ierr = PetscFree(self->ts_s);  CHKERRQ(ierr);
  ierr = PetscFree(self->ts_r);  CHKERRQ(ierr);  
  ierr = PetscFree(self->xs); CHKERRQ(ierr);
  ierr = PetscFree(self->xs_s); CHKERRQ(ierr);
  ierr = PetscFree(self->ws); CHKERRQ(ierr);
  ierr = PetscFree(self->qrs); CHKERRQ(ierr);
  ierr = PetscFree(self->Rrs); CHKERRQ(ierr);
  
  ierr = PetscFree(self->vals); CHKERRQ(ierr);
  ierr = PetscFree(self->derivs); CHKERRQ(ierr);
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BSSView(BSS self, PetscViewer v) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

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
    PetscViewerASCIIPrintf(v, "BSS object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "order: %d\n", self->order);  
    PetscViewerASCIIPrintf(v, "num_ele: %d\n", self->num_ele);
    PetscViewerASCIIPrintf(v, "num_basis: %d\n", self->num_basis);
    BPSView(self->bps, v);  
    CScalingView(self->c_scaling, v);
    if(format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      PetscViewerASCIIPrintf(v, "b_idx_list: ");
      for(int i = 0; i < self->num_basis; i++)
	PetscViewerASCIIPrintf(v, " %d ", self->b_idx_list[i]);	
    }
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}
PetscErrorCode BSSCheck(BSS self) {
  if(self == NULL)
    SETERRQ(self->comm, 1, "object is null");
  
  if(!self->setup)
    SETERRQ(self->comm, 1, "BSS object is not setup. Call BSSSetKnots and BSSSetUp first");

  return 0;
}

PetscErrorCode BSSSetKnots(BSS self, int order, BPS bps) {

  PetscErrorCode ierr;
  

  self->order = order;
  self->bps = bps;
  ierr = BPSGetNumEle(bps, &self->num_ele); CHKERRQ(ierr);

  PetscReal *zs; PetscInt num_zs;
  BPSGetZs(self->bps, &zs, &num_zs);
  
  self->num_basis = NumBSpline(self->order, num_zs-1);
  self->num_ts = num_zs+2*self->order-2;

  PetscMalloc1(self->num_basis, &self->b_idx_list);

  PetscMalloc1(self->num_ts, &self->ts_r);
  PetscMalloc1(self->num_ts, &self->ts_s);

  int n_xs = self->num_ele * self->order;
  PetscMalloc1(n_xs, &self->xs); 
  PetscMalloc1(n_xs, &self->xs_s); 
  PetscMalloc1(n_xs, &self->ws);
  PetscMalloc1(n_xs, &self->qrs); 
  PetscMalloc1(n_xs, &self->Rrs);

  int num = n_xs*self->num_basis;
  PetscMalloc1(num, &self->vals); 
  PetscMalloc1(num, &self->derivs);


  PetscFree(zs);
  self->set_knots = PETSC_TRUE;
  return 0;
}
PetscErrorCode BSSSetCScaling(BSS self, CScaling c_scaling) {

  CScaling s;
  if(c_scaling == NULL) {
    CScalingCreate(self->comm, &s); CScalingSetNone(s);
  } else {
    s = c_scaling;
  }
  self->c_scaling = s;
  return 0;
}
PetscErrorCode BSSSetUp(BSS self) {

  PetscErrorCode ierr;

  // check knots are setted
  if(!self->set_knots)
    SETERRQ(self->comm, 1, "Knots information is not setted. Call BSSSetKnots first.");

  //  set default scaler
  if(self->c_scaling == NULL) 
    BSSSetCScaling(self, NULL);
  
  // copy ts_r and  ts_s
  PetscReal *zs; PetscInt num_zs;
  BPSGetZs(self->bps, &zs, &num_zs);

  for(int i = 0; i < self->order-1; i++) {
    self->ts_r[i] = zs[0];
    self->ts_r[self->order-1+num_zs+i] = zs[num_zs-1];
  }
  for(int i = 0; i < num_zs; i++) 
    self->ts_r[i+self->order-1] = zs[i];
  PetscFree(zs);
  CScalingCalc(self->c_scaling, self->ts_r, self->num_ts, NULL, self->ts_s);
  //  for(int i = 0; i < self->num_ts; i++)
  //    self->ts_s[i] = self->ts_r[i];

  // index of basis
  for(int ib = 0; ib < self->num_basis; ib++)
    self->b_idx_list[ib] = ib + 1;

  // each element
  for(int ie = 0; ie < self->num_ele; ie++) {
    PetscScalar a, b; a = self->bps->zs[ie]; b = self->bps->zs[ie+1];
    for(int iq = 0; iq < self->order; iq++) {
      PetscReal x_r[1], w_r[1];
      PetscScalar x_c[1];
      int ix = ie*self->order+iq;
      PetscScalar quad_x, quad_w;
      LegGauss(self->order, iq, &quad_x, &quad_w);
      x_r[0] = creal(quad_x); w_r[0] = creal(quad_w);
      x_r[0] = (b+a)/2.0  + (b-a)/2.0 * x_r[0]; w_r[0] = (b-a)/2.0*w_r[0];
      CScalingCalc(self->c_scaling, x_r, 1, NULL, x_c);
      self->xs[ix] = x_r[0]; 
      self->ws[ix] = w_r[0];
      self->xs_s[ix] = x_c[0];
            
      for(int ib = 0; ib < self->num_basis; ib++) {
	PetscScalar y;
	PetscScalar dy;
	PetscReal *ts_r = self->ts_r;
	PetscScalar *ts_s = self->ts_s;
	int idx = self->b_idx_list[ib];
	CalcBSpline(     self->order, ts_r, ts_s, idx, x_r[0], x_c[0], &y,  NULL);
	CalcDerivBSpline(self->order, ts_r, ts_s, idx, x_r[0], x_c[0], &dy, NULL);
	int iy = ib*(self->num_ele*self->order) + ie*self->order + iq;
	self->vals[iy] = y;  
	self->derivs[iy] = dy;
      }
    }
  }

  int n_xs = self->num_ele * self->order;
  ierr = CScalingCalc(self->c_scaling, self->xs, n_xs,
		      self->qrs, self->Rrs); CHKERRQ(ierr);
  self->setup = PETSC_TRUE;
  
  return 0;
}
PetscErrorCode BSSSetFromOptions(BSS self) {

  PetscErrorCode ierr;
  PetscBool find;

  PetscInt order = 2;
  BPS bps;        BPSCreate(self->comm, &bps);
  CScaling scaler;  CScalingCreate(self->comm, &scaler);
  
  ierr = PetscOptionsGetInt(NULL, NULL, "-bss_order", &order, &find); CHKERRQ(ierr);
  ierr = BPSSetFromOptions(bps); CHKERRQ(ierr);
  ierr = CScalingSetFromOptions(scaler); CHKERRQ(ierr);
  ierr = BSSSetKnots(self, order, bps); CHKERRQ(ierr);
  ierr = BSSSetCScaling(self, scaler); CHKERRQ(ierr);
  ierr = BSSSetUp(self); CHKERRQ(ierr);

  return 0;
 }

PetscErrorCode BSSPsiOne(BSS self, Vec c, PetscReal x, PetscScalar *y) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  int n; BSSGetSize(self, &n);
  Vec us; BSSCreateR1Vec(self, &us);
  for(int i = 0; i < n; i++) {
    PetscScalar u;
    ierr = BSSBasisPsi(self, i, x, &u); CHKERRQ(ierr);
    VecSetValue(us, i, u, INSERT_VALUES);
  }
  VecAssemblyBegin(us); VecAssemblyEnd(us);

  PetscScalar yy;
  VecTDot(us, c, &yy);
  *y = yy;
  VecDestroy(&us);
  return 0;
}
PetscErrorCode BSSPsi(BSS self, Vec c, Vec xs, Vec ys) {
  
  // yj = sum(i) ci ui(xj)

  //  PetscErrorCode ierr;
  //  ierr = BSSChcek(self); CHKERRQ(ierr);
  PetscErrorCode ierr;  
  int nc; VecGetSize(c, &nc);
  int nb; BSSGetSize(self, &nb);
  int nx; VecGetSize(xs, &nx);
  int ny; VecGetSize(ys, &ny);
  if(nc != nb) {
    SETERRQ(self->comm, 1, "size of c must be same as basis size");
  }
  if(nx != ny) {
    SETERRQ(self->comm, 1, "xs and ys must be same size");
  }

  Mat f_jx_ib;
  ierr = MatCreate(self->comm, &f_jx_ib); CHKERRQ(ierr);
  ierr = MatSetSizes(f_jx_ib, PETSC_DECIDE, PETSC_DECIDE, nx, nb); CHKERRQ(ierr);
  ierr = MatSetUp(f_jx_ib);               CHKERRQ(ierr);

  PetscScalar *x_ptr;
  VecGetArray(xs, &x_ptr);

  PetscReal *x_r_ptr; PetscMalloc1(nx, &x_r_ptr);
  for(int ix = 0; ix < nx; ix++)
    x_r_ptr[ix] = creal(x_ptr[ix]);

  PetscScalar *x_c_ptr; PetscMalloc1(nx, &x_c_ptr);
  ierr = CScalingCalc(self->c_scaling, x_r_ptr, nx, NULL, x_c_ptr); CHKERRQ(ierr);

  for(int jx = 0; jx < nx; jx++) {
    for(int ib = 0; ib < nb; ib++) {
      PetscScalar y;
      PetscBool zeroq;
      ierr = CalcBSpline(self->order, self->ts_r, self->ts_s, self->b_idx_list[ib],
			 x_r_ptr[jx], x_c_ptr[jx], &y, &zeroq); CHKERRQ(ierr);
      if(!zeroq) {
	ierr = MatSetValue(f_jx_ib, jx, ib, y, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }  
  MatAssemblyBegin(f_jx_ib, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(f_jx_ib, MAT_FINAL_ASSEMBLY);

  ierr = MatMult(f_jx_ib, c, ys); CHKERRQ(ierr);

  MatDestroy(&f_jx_ib);
  VecRestoreArray(xs, &x_ptr);
  PetscFree(x_r_ptr);
  PetscFree(x_c_ptr);
  return 0;
}
PetscErrorCode BSSDerivPsiOne(BSS self, Vec c, PetscReal x, PetscScalar *y) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  int n;
  ierr = BSSGetSize(self, &n); CHKERRQ(ierr);
  Vec us;
  ierr = BSSCreateR1Vec(self, &us); CHKERRQ(ierr);

  for(int i = 0; i < n; i++) {
    PetscScalar u;
    ierr = BSSDerivBasisPsi(self, i, x, &u); CHKERRQ(ierr);
    ierr = VecSetValue(us, i, u, INSERT_VALUES); CHKERRQ(ierr);
  }
  VecAssemblyBegin(us); VecAssemblyEnd(us);

  PetscScalar yy;
  ierr = VecTDot(us, c, &yy); CHKERRQ(ierr);
  *y = yy;
  VecDestroy(&us);
  return 0;
}
PetscErrorCode BSSDerivPsi(BSS self, Vec c, Vec xs, Vec ys) {
  PetscErrorCode ierr;  
  int nc; VecGetSize(c, &nc);
  int nb; BSSGetSize(self, &nb);
  int nx; VecGetSize(xs, &nx);
  int ny; VecGetSize(ys, &ny);
  if(nc != nb) {
    SETERRQ(self->comm, 1, "size of c must be same as basis size");
  }
  if(nx != ny) {
    SETERRQ(self->comm, 1, "xs and ys must be same size");
  }

  Mat f_jx_ib;
  ierr = MatCreate(self->comm, &f_jx_ib); CHKERRQ(ierr);
  ierr = MatSetSizes(f_jx_ib, PETSC_DECIDE, PETSC_DECIDE, nx, nb); CHKERRQ(ierr);
  ierr = MatSetUp(f_jx_ib);               CHKERRQ(ierr);

  PetscScalar *x_ptr;
  VecGetArray(xs, &x_ptr);

  PetscReal *x_r_ptr; PetscMalloc1(nx, &x_r_ptr);
  for(int ix = 0; ix < nx; ix++)
    x_r_ptr[ix] = creal(x_ptr[ix]);

  PetscScalar *x_c_ptr; PetscMalloc1(nx, &x_c_ptr);
  ierr = CScalingCalc(self->c_scaling, x_r_ptr, nx, NULL, x_c_ptr); CHKERRQ(ierr);

  for(int jx = 0; jx < nx; jx++) {
    for(int ib = 0; ib < nb; ib++) {
      PetscScalar y;
      PetscBool zeroq;
      ierr = CalcDerivBSpline(self->order,
			      self->ts_r, self->ts_s, self->b_idx_list[ib],
			      x_r_ptr[jx], x_c_ptr[jx], &y, &zeroq); CHKERRQ(ierr);
      if(!zeroq) {
	ierr = MatSetValue(f_jx_ib, jx, ib, y, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }  
  MatAssemblyBegin(f_jx_ib, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(f_jx_ib, MAT_FINAL_ASSEMBLY);

  ierr = MatMult(f_jx_ib, c, ys); CHKERRQ(ierr);

  MatDestroy(&f_jx_ib);
  VecRestoreArray(xs, &x_ptr);
  PetscFree(x_r_ptr);
  PetscFree(x_c_ptr);
  return 0;  
}
PetscErrorCode BSSBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  PetscScalar z;

  PetscReal xs[1] = {x};
  PetscScalar Rxs[1];
  CScalingCalc(self->c_scaling, xs, 1, NULL, Rxs);
  
  CalcBSpline(self->order, 
	      self->ts_r, 
	      self->ts_s, 
	      self->b_idx_list[i], 
	      xs[0], 
	      Rxs[0],
	      &z,
	      NULL);
  *y = z;
  return 0;
}
PetscErrorCode BSSDerivBasisPsi(BSS self, int i, PetscReal x, PetscScalar *y) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  PetscReal xs[1] = {x};
  PetscScalar Rxs[1];
  CScalingCalc(self->c_scaling, xs, 1, NULL, Rxs);

  PetscScalar z;
  CalcDerivBSpline(self->order, 
		   self->ts_r, 
		   self->ts_s, 
		   self->b_idx_list[i], 
		   xs[0], 
		   Rxs[0],
		   &z,
		   NULL);
  *y = z;
  return 0;
}
PetscErrorCode BSSGetSize(BSS self, int *n) {
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  *n = self->num_basis;
  return 0;
}
 
// ---- Matrix -----
PetscErrorCode BSSCreateR1Mat(BSS self, Mat *M) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  int nb = self->num_basis;
  ierr = MatCreate(self->comm, M);  CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nb, nb);  CHKERRQ(ierr);
  ierr = MatSetUp(*M);  CHKERRQ(ierr);
  return 0;
}
PetscErrorCode BSSCreateR2Mat(BSS self, Mat *M) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  int nb = self->num_basis;
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, nb*nb, nb*nb);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode BSSCreateR1Vec(BSS self, Vec *v) {

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  int nb = self->num_basis;
  VecCreate(self->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, nb);
  VecSetUp(*v);
  return 0;
}

PetscErrorCode BSSSR1Mat(BSS self, Mat M) {
  int i, j, k;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = self->order;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  InsertMode mode = INSERT_VALUES;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(self->order, self->b_idx_list[i], self->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++) 
	  v += self->vals[k+i*(ne*nq)] * self->vals[k+j*(ne*nq)] * self->ws[k] * self->qrs[k];
	//v += self->vals[k+i*(ne*nq)] * self->vals[k+j*(ne*nq)] * self->ws[k];
	ierr = MatSetValue(M, i, j, v, mode); CHKERRQ(ierr);
      }
    }

  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSR2invR1Mat(BSS self, Mat M) {

  int i, j, k;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = self->order;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(self->order, self->b_idx_list[i], self->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++) {
	  PetscScalar x = self->Rrs[k];
	  PetscReal w = self->ws[k];
	  PetscScalar q = self->qrs[k];
	  PetscScalar fi = self->vals[k+i*ne*nq];
	  PetscScalar fj = self->vals[k+j*ne*nq];
	  v += fi*fj*w*q/(x*x);
	}
	ierr = MatSetValue(M, i, j, v, mode); CHKERRQ(ierr);
      }
    }

  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSD2R1Mat(BSS self, Mat M) {
  int i, j, k;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = self->order;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  for(i = 0; i < nb; i++)
    for(j = 0; j < nb; j++) {
      if(HasNon0Value(self->order, self->b_idx_list[i], self->b_idx_list[j])) {
	PetscScalar v = 0.0;
	for(k = 0; k < ne*nq; k++)  {
	  /*
	  v += (self->derivs[k+i*(ne*nq)] *
		self->derivs[k+j*(ne*nq)] *
		self->ws[k] / (self->qrs[k]*self->qrs[k]));
	  */
	  /*
	  v += (self->derivs[k+i*(ne*nq)] *
		self->derivs[k+j*(ne*nq)] *
		self->ws[k]);
	  */
	  v += (self->derivs[k+i*(ne*nq)] *
		self->derivs[k+j*(ne*nq)] *
		self->ws[k] * self->qrs[k]);
	  /*
	  v += (self->derivs[k+i*(ne*nq)] *
		self->derivs[k+j*(ne*nq)] *
		self->ws[k] / self->qrs[k]);
	  */
	  //v += self->derivs[k+i*(ne*nq)] * self->derivs[k+j*(ne*nq)] * self->ws[k] / self->qrs[k];
	  //	  v += self->derivs[k+i*(ne*nq)] * self->derivs[k+j*(ne*nq)] * self->ws[k];
	}
	
	ierr = MatSetValue(M, i, j, -v, mode); CHKERRQ(ierr);
      }
    }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSENR1Mat(BSS self, int q, PetscReal a, Mat M) {
  int order = self->order;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = order*ne;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  PetscScalar *vs; 
  PetscMalloc1(nq, &vs);
  for(int k = 0; k < nq; k++) {
    PetscScalar g = self->xs[k] > a ? self->Rrs[k] : a;
    PetscScalar s = self->xs[k] < a ? self->Rrs[k] : a;
    //vs[k] = pow(s/g, q)/g * self->ws[k];
    vs[k] = pow(s/g, q)/g * self->ws[k] * self->qrs[k];
  }

  for(int i = 0; i < nb; i++) {
    int j0 = i-order+1; j0 = j0 < 0 ? 0: j0;
    for(int j = j0; j <= i; j++) {
      int k0, k1;  
      Non0QuadIndex(j, i, order, nq, &k0, &k1);
      PetscScalar v = 0.0;
      for(int k = k0; k < k1; k++) 
	v += self->vals[k+i*nq] * self->vals[k+j*nq] * vs[k];
      ierr = MatSetValue(M, i, j, v, mode); CHKERRQ(ierr);
      if(i!=j)
	ierr = MatSetValue(M, j, i, v, mode); CHKERRQ(ierr);
    }
  }

  PetscFree(vs);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSPotR1Mat(BSS self, PF pot, Mat M) {
  int order = self->order;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = order*ne;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  PetscScalar *vs; PetscMalloc1(nq, &vs);
  //ierr = PFApply(pot, nq, self->xs_s, vs); CHKERRQ(ierr);
  ierr = PFApply(pot, nq, self->Rrs, vs); CHKERRQ(ierr);
  for(int k = 0; k < nq; k++) {
    vs[k] *= self->ws[k] * self->qrs[k];
  }

  for(int i = 0; i < nb; i++) {
    int j0 = i-order+1; j0 = j0 < 0 ? 0: j0;
    for(int j = j0; j <= i; j++) {
      int k0, k1;  
      Non0QuadIndex(j, i, order, nq, &k0, &k1);
      PetscScalar v = 0.0;
      for(int k = k0; k < k1; k++) 
	v += self->vals[k+i*nq] * self->vals[k+j*nq] * vs[k];
      ierr = MatSetValue(M, i, j, v, mode); CHKERRQ(ierr);
      if(i!=j)
	ierr = MatSetValue(M, j, i, v, mode); CHKERRQ(ierr);
    }
  }

  PetscFree(vs);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;  
}
PetscErrorCode BSSPotR1Vec(BSS self, PF pot, Vec V) {
  int order = self->order;
  int nb = self->num_basis;
  int ne = self->num_ele;
  int nq = order*ne;

  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);

  PetscScalar *vs; PetscMalloc1(nq, &vs);  
  ierr = PFApply(pot, nq, self->Rrs, vs); CHKERRQ(ierr);
  for(int i = 0; i < nb; i++) {
    int k0, k1; Non0QuadIndex(i, i, order, nq, &k0, &k1);
    PetscScalar v = 0.0;
    for(int k = k0; k < k1; k++) {
      //v += self->vals[k+i*nq] * self->ws[k] * vs[k];
      v += self->vals[k+i*nq] * self->ws[k] * self->qrs[k] * vs[k];
    }
    ierr = VecSetValue(V, i, v, INSERT_VALUES); CHKERRQ(ierr);
  }

  PetscFree(vs);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);

  return 0;
}
PetscErrorCode BSSOpMat(BSS self, Op op, Mat M) {

  SETERRQ(self->comm, 1, "not implemented");

  /*
  PetscErrorCode ierr;
  if(OpIsType(op, OpD2)) {
    ierr = BSSD2R1Mat(self, M); CHKERRQ(ierr);
  } else if(OpIsType(op, OpPF)) {
    PF pf;
    ierr = OpGetPF(op, &pf);
    ierr = BSSPotR1Mat(self, pf, M); CHKERRQ(ierr);
  } else if(OpIsType(op, OpPartialVee)) {
  SETERRQ(self->comm, 1, "unsupported tpye");
  } else 
    SETERRQ(self->comm, 1, "unsupported tpye");

  return 0;
  */
}
PetscErrorCode BSSEER2Mat(BSS self, int q, Mat M) {

  int k = self->order;
  int nb = self->num_basis;
  int nq = self->order * self->num_ele;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  PetscScalar *sg_ij; 
  ierr = PetscMalloc1(nq*nq, &sg_ij); CHKERRQ(ierr);
  for(int i = 0; i < nq; i++)
    for(int j = 0; j < nq; j++) {
      PetscScalar g = self->xs[i]>self->xs[j] ? self->Rrs[i] : self->Rrs[j];
      PetscScalar s = self->xs[i]>self->xs[j] ? self->Rrs[j] : self->Rrs[i];
      sg_ij[j+nq*i] = pow(s/g, q)/g;      
    }

  PetscScalar *wvv; ierr = PetscMalloc1(nb*nb*nq, &wvv); CHKERRQ(ierr);
  int *i0s; ierr = PetscMalloc1(nb*nb, &i0s); CHKERRQ(ierr);
  int *i1s; ierr = PetscMalloc1(nb*nb, &i1s); CHKERRQ(ierr);
  for(int a = 0; a < nb; a++){
    int c0 = a-k+1; c0 = c0<0?0:c0;
    int c1 = a+k  ; c1 = c1>nb?nb:c1;
    for(int c = c0; c < c1; c++) {
      int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
      i0s[a*nb+c] = i0; i1s[a*nb+c] = i1;
      for(int n = 0; n < nq; n++) 
	wvv[a*nb*nq+c*nq+n]= self->ws[n] * self->vals[a*nq+n] * self->vals[c*nq+n];
    }
  }

  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1; c0 = c0<0?0:c0;
    for(int c = c0; c <= a; c++) {
      for(int b = 0; b < nb; b++) {
	int d0 = b-k+1; d0 = d0<0?0:d0;
	for(int d = d0; d <= b; d++) {
	  PetscScalar v = 0.0;
	  for(int i = i0s[a*nb+c]; i < i1s[a*nb+c]; i++) {
	    PetscScalar aci = wvv[a*nb*nq + c*nq + i];
	    for(int j = i0s[b*nb+d]; j < i1s[b*nb+d]; j++)
	      v += aci * sg_ij[i*nq+j] * wvv[b*nb*nq + d*nq + j];
	  }
	  
	  if(a > c && b > d) {
	    ierr = MatSetValue(M, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(M, a*nb+d, c*nb+b, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(M, c*nb+b, a*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(M, c*nb+d, a*nb+b, v, mode); CHKERRQ(ierr);
	  } else if(a > c && b == d) {
	    ierr = MatSetValue(M, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(M, c*nb+b, a*nb+d, v, mode); CHKERRQ(ierr);
	  } else if(a == c && b > d) {
	    ierr = MatSetValue(M, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	    ierr = MatSetValue(M, a*nb+d, c*nb+b, v, mode); CHKERRQ(ierr);   
	  } else if(a == c && b == d) {
	    ierr = MatSetValue(M, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	  }
	}
      }
    }
  }

  PetscFree(sg_ij); PetscFree(wvv); PetscFree(i0s); PetscFree(i1s);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode BSSEER2Mat_ver1(BSS self, int q, Mat M) {

  int nb = self->num_basis;
  int nq = self->order * self->num_ele;
  PetscErrorCode ierr;
  ierr = BSSCheck(self); CHKERRQ(ierr);
  InsertMode mode = INSERT_VALUES;

  double *sg_ij;
  sg_ij = (double*)malloc(sizeof(double)*nq*nq);
  for(int i = 0; i < nq; i++)
    for(int j = 0; j < nq; j++) {
      double v; PartialCoulomb(q, self->xs[i], self->xs[j], &v);
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
	      v += self->ws[n] * self->vals[a*nq+m] * self->vals[c*nq+m] 
		* self->ws[m] * self->vals[b*nq+n] * self->vals[d*nq+n] 
		* sg_ij[n*nq+m];
	    }
	  }
	  ierr = MatSetValue(M, a*nb+b, c*nb+d, v, mode); CHKERRQ(ierr);
	}
      }
    }
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
