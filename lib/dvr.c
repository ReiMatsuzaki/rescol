#include "../include/math.h"
#include "../include/mat.h"
#include "../include/synthesize.h"
#include "../include/dvr.h"

// ------- Lapack -------
int dgetrf_(long*, long*, PetscReal*,   long*, long*, long*);
int zgetrf_(long*, long*, PetscScalar*, long*, long*, long*);
int dgetri_(long*, double*, long*, long*, double*, long*, long*);
int zgetri_(long*, PetscScalar*, long*, long*, PetscScalar*, long*, long*);

// ------- external functions ---------
PetscErrorCode MatCreateTransformMat(PetscScalar *ws_c, int nq, int ne, 
				     MPI_Comm comm, Mat *A) {
  MatCreate(comm, A);
  MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, nq*ne, (nq-2)*ne + ne-1);
  MatSetUp(*A);

  int ii = 1;
  int jj = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 1; i_q < nq-1; i_q++) {
      MatSetValue(*A, ii, jj, 1.0/csqrt(ws_c[i_ele*nq+i_q]), INSERT_VALUES);
      ii++; jj++;
    }
    if(i_ele != ne-1) {
      PetscScalar d = 1.0/csqrt(ws_c[i_ele*nq+nq-1] + ws_c[(1+i_ele)*nq]);
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
PetscErrorCode ValueLS(BPS bps, PetscScalar *xs_c, PetscInt nq, 
		       PetscInt i, PetscInt m, PetscReal x,
		       PetscScalar x_c, PetscScalar *v, PetscBool *zeroq) {
  /*
    Gives value of Lobatto shape function f_{i,m} at quadrature point x_{i,mp}.

    bps : break points set 
    xs_c: quadrature points (Scalar[nq])
    nq : number of quadrature in each element
    i : index of element which represent LS function
    m : index of quadrature point which represent  Lobatto Shape function
    x : target location
  */

  // if x is out of ith element.
  PetscBool is_in_q;
  BPSInElementQ(bps, i, x, &is_in_q);
  if(!is_in_q) {
    *zeroq = PETSC_TRUE;
    *v = 0.0;
    return 0;
  }

  *zeroq = PETSC_FALSE;
  // compute Lobatto Shape function
  PetscScalar cumsum = 1.0;
  for(int iq = 0; iq < nq; iq++) {
    if(iq != m) {
      //printf("%f\n", creal(x_c-xs_c[i*nq + iq])/(xs_c[i*nq+m]-xs_c[i*nq + iq]));
      cumsum *= (x_c-xs_c[i*nq + iq])/(xs_c[i*nq+m]-xs_c[i*nq + iq]);
    }
  }
  *v = cumsum;
  return 0;
}
PetscErrorCode ValueDerivLS(BPS bps, PetscScalar *xs_c, PetscInt nq, 
			    PetscInt i, PetscInt m, PetscReal x,
			    PetscScalar x_c, PetscScalar *v, PetscBool *zeroq) {
  /*
    Gives value of derivative of Lobatto shape function f_{i,m} at 
    quadrature point x_{i,mp}.

    xs_c: quadrature points
    zs  : knot points (Scalar[ne+1])
    nq : number of quadrature in each element
    ne : number of element
    i : index of element which represent Lobatto Shape function
    m : index of quadrature point which represent  Lobatto Shape function
    x : target location
  */

  // if x is out of ith element.
  PetscBool is_in_q; BPSInElementQ(bps, i, x, &is_in_q);
  if(!is_in_q) {
    *v = 0.0;
    *zeroq = PETSC_TRUE;
    return 0;
  }
  *zeroq = PETSC_FALSE;
  // compute Lobatto Shape function
  PetscScalar cumsum = 0.0;
  for(int j = 0; j < nq; j++) {
    if(j != m) {
      PetscScalar cumprod = 1.0/(xs_c[i*nq+m]-xs_c[i*nq + j]);
      for(int iq = 0; iq < nq; iq++) {
	if(iq != m && iq != j) {
	  cumprod *= (x_c-xs_c[i*nq+iq])/(xs_c[i*nq+m]-xs_c[i*nq+iq]);
	}
      }
      cumsum += cumprod;
    }
    
  }
  *v = cumsum;
  return 0;
}
PetscErrorCode DerivLS(PetscScalar *xs_c, PetscScalar *ws_c, 
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
      *v = +1.0/(2.0*ws_c[i*nq+m]);
    else if(m == 0) 
      *v = -1.0/(2.0*ws_c[i*nq+m]);
    else
      *v = 0.0;
  } else {
    *v = 1.0/(xs_c[i*nq+m]-xs_c[i*nq+mp]);
    for (int k = 0; k < nq; k++) 
      if(k != m && k != mp) 
	*v *= (xs_c[i*nq+mp]-xs_c[i*nq+k])/(xs_c[i*nq+m]-xs_c[i*nq+k]);
  }
  return 0;
}

// ---- Basic Methods ----
PetscErrorCode DVRCreate(MPI_Comm comm, DVR *p_self) {
  DVR self;
  PetscMalloc1(1, &self);

  self->comm = comm;
  self->nq = -1;
  self->bps = NULL;
  CScalingCreate(comm, &self->cscaling);
  CScalingSetNone(self->cscaling);

  self->num_basis = -1;
  self->xs = NULL;
  self->xs_basis = NULL;
  self->xs_c = NULL;
  self->ws_c = NULL;
  self->xs_c = NULL;
  self->ws_c = NULL;
  self->xs_basis_c = NULL;
  self->ws_basis_c = NULL;
  self->vs_basis = NULL;

  self->D2_R1LSMat = NULL;
  self->R2_R1LSMat = NULL;
  self->T = NULL;
  self->TT = NULL;
  self->T2 = NULL;
  self->T2T = NULL;

  *p_self = self;
  return 0;
}
PetscErrorCode DVRCopy(DVR self, DVR other) {

  PetscErrorCode ierr;

  BPS bps;
  ierr = BPSCreate(self->comm, &bps); CHKERRQ(ierr);
  ierr = BPSCopy(self->bps, bps); CHKERRQ(ierr);
  ierr = DVRSetKnots(other, self->nq, bps); CHKERRQ(ierr);

  CScaling cscaling;
  ierr = CScalingCreate(self->comm, &cscaling); CHKERRQ(ierr);
  ierr = CScalingCopy(self->cscaling, cscaling); CHKERRQ(ierr);
  ierr = DVRSetCScaling(other, cscaling); CHKERRQ(ierr);

  ierr = DVRSetUp(other); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode DVRDestroy(DVR *p_self) {
  PetscErrorCode ierr;
  DVR self = *p_self;
  ierr = BPSDestroy(&self->bps);  CHKERRQ(ierr);    
  ierr = CScalingDestroy(&self->cscaling); CHKERRQ(ierr);
  ierr = PetscFree(self->xs); CHKERRQ(ierr);    
  ierr = PetscFree(self->xs_basis); CHKERRQ(ierr);    
  ierr = PetscFree(self->xs_c); CHKERRQ(ierr);    
  ierr = PetscFree(self->ws_c); CHKERRQ(ierr);    
  ierr = PetscFree(self->xs_basis_c); CHKERRQ(ierr);
  ierr = PetscFree(self->ws_basis_c); CHKERRQ(ierr);
  ierr = PetscFree(self->vs_basis); CHKERRQ(ierr);

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
    PetscViewerASCIIPrintf(v, "DVR object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "nq: %d\n", self->nq);
    CScalingView(self->cscaling, v);
    
    PetscViewerASCIIPrintf(v, "num_basis: %d\n", self->num_basis);
    BPSView(self->bps, v);  
    int n_basis; DVRGetSize(self, &n_basis);

    /*
    PetscViewerASCIIPrintf(v, "xs_basis | xs_basis_c | ws_basis_c:\n");
    for(int i = 0; i < n_basis; i++) {
      PetscViewerASCIIPrintf(v, "%f | %f, %f | %f, %f \n",
			     creal(self->xs_basis[i]),
			     creal(self->xs_basis_c[i]),
			     cimag(self->xs_basis_c[i]),
			     creal(self->ws_basis_c[i]),
			     cimag(self->ws_basis_c[i]));
    }
    int ne; BPSGetNumEle(self->bps, &ne);
    int nq = self->nq;
    PetscViewerASCIIPrintf(v, "xs:\n");
    for(int i = 0; i < ne*nq; i++) {
      PetscViewerASCIIPrintf(v, "%f | %f, %f | %f | %f, %f \n",
			     creal(self->xs[i]),
			     creal(self->xs_c[i]),
			     cimag(self->xs_c[i]),
			     self->ws[i],
			     creal(self->ws_c[i]),
			     cimag(self->ws_c[i])
			     );
    }
    */
    if(format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      PetscViewerASCIIPrintf(v, "xs: ");
      for(int i = 0; i < self->nq; i++)
	PetscViewerASCIIPrintf(v, " %f ", self->xs[i]);	
    }
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}

// ---- Accessor ----
PetscErrorCode DVRSetKnots(DVR self, int nq, BPS bps) {

  if(bps == NULL) {
    SETERRQ(self->comm, 1, "bps is numm");
  }

  self->nq = nq;
  self->bps = bps;

  return 0;
}
PetscErrorCode DVRSetCScaling(DVR self, CScaling _cscaling) {
  CScalingDestroy(&self->cscaling);
  self->cscaling = _cscaling;
  return 0;
}
PetscErrorCode DVRSetUp(DVR self) {

  PetscErrorCode ierr;

  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int nx = nq * ne;
  ierr = PetscMalloc1(nx, &self->xs); CHKERRQ(ierr);
  ierr = PetscMalloc1(nx, &self->xs_c); CHKERRQ(ierr);
  ierr = PetscMalloc1(nx, &self->ws_c); CHKERRQ(ierr);
  ierr = PetscMalloc1(nx, &self->ws_basis_c); CHKERRQ(ierr);
  self->num_basis = NumDVRBasis(self->nq, ne);  
  ierr = PetscMalloc1(self->num_basis, &self->xs_basis); CHKERRQ(ierr);
  ierr = PetscMalloc1(self->num_basis, &self->xs_basis_c); CHKERRQ(ierr);
  ierr = PetscMalloc1(self->num_basis, &self->vs_basis); CHKERRQ(ierr);
  
  for(int i = 0; i < self->num_basis; i++) {
    self->vs_basis[i] = 777.0;
  }

  PetscReal R0;    CScalingGetRadius(self->cscaling, &R0);
  PetscReal theta; CScalingGetTheta(self->cscaling, &theta);
  PetscBool use_cscaling; CScalingQ(self->cscaling, &use_cscaling);
  
  double eps = 0.0000001;
  PetscScalar scale_factor = cexp(I*theta*M_PI/180.0);
  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    PetscReal a, b;
    BPSGetEdge(self->bps, i_ele, &a, &b);
    PetscBool do_scale = (use_cscaling &&
			  R0-eps<a     &&
			  R0-eps<b);
    for(int i_q = 0; i_q < nq; i_q++) {
      PetscScalar x, w;
      LobGauss(self->nq, i_q, &x, &w);
      x = (b+a)/2.0 + (b-a)/2.0*x;
      w = (b-a)/2.0*w;
      self->xs[idx] = x; //self->ws[idx] = w;
      if(do_scale) {
	CScalingCalcOne(self->cscaling, x, &scale_factor,
			&self->xs_c[idx]);
	self->ws_c[idx] = w*scale_factor;
      } else {
	self->xs_c[idx] = x;
	self->ws_c[idx] = w;
      }
      idx++;
    }
  }
  idx = 0;
  self->vs_basis[77] = 77.0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 1; i_q < nq-1; i_q++) {
      self->xs_basis[idx] = self->xs[i_ele*nq+i_q];      
      self->xs_basis_c[idx] = self->xs_c[i_ele*nq+i_q];  
      // modified @ 2016/6/21
      //self->ws_basis_c[idx] = self->ws_c[(i_ele+1)*nq];
      self->ws_basis_c[idx] = self->ws_c[i_ele*nq+i_q];
      self->vs_basis[idx] = 1.0/sqrt(self->ws_basis_c[idx]);
      idx++;      
    }
    if(i_ele != ne-1) {
      self->xs_basis[idx] = self->xs[(i_ele+1)*nq];
      self->xs_basis_c[idx] = self->xs_c[(i_ele+1)*nq];
      self->ws_basis_c[idx] = self->ws_c[(i_ele+1)*nq];      
      self->vs_basis[idx]
	= 1.0/sqrt(self->ws_c[(i_ele+1)*nq] + self->ws_c[i_ele*nq+(nq-1)]);
      idx++;
    }
  }
  self->D2_R1LSMat = NULL;
  self->R2_R1LSMat = NULL;
  ierr = MatCreateTransformMat(self->ws_c, nq, ne, self->comm, &self->T); 
  CHKERRQ(ierr);
  ierr = MatTranspose(self->T, MAT_INITIAL_MATRIX, &self->TT); CHKERRQ(ierr);
  self->T2 = NULL;
  self->T2T = NULL;

  return 0;
}
PetscErrorCode DVRSetFromOptions(DVR self) {

  PetscBool find;
  PetscInt nq = 2;  
  PetscErrorCode ierr;

  BPS bps;
  ierr = BPSCreate(self->comm, &bps); CHKERRQ(ierr);
  ierr = BPSSetFromOptions(bps); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-dvr_nq", &nq, &find); CHKERRQ(ierr);
  ierr = DVRSetKnots(self, nq, bps); CHKERRQ(ierr);

  CScaling cscaling;
  ierr = CScalingCreate(self->comm, &cscaling); CHKERRQ(ierr);
  ierr = CScalingSetFromOptions(cscaling); CHKERRQ(ierr);
  ierr = DVRSetCScaling(self, cscaling);   CHKERRQ(ierr);

  ierr = DVRSetUp(self); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode DVRBasisPsi(DVR self, int i, PetscScalar x, PetscScalar *y) {
  SETERRQ(self->comm, 1, "not implemented yet");
  *y = i + x;
  return 0;
}
PetscErrorCode DVRGetSize(DVR self, int *n) {
  *n = self->num_basis;
  return 0;
}
PetscErrorCode DVRGetLSSize(DVR self, int *n) {
  int ne; BPSGetNumEle(self->bps, &ne);
  *n = ne * self->nq;
  return 0;
}

// ----  Calculation ----
PetscErrorCode DVRPsiOne(DVR self, Vec cs, PetscReal x, PetscScalar *y) {

  PetscScalar xc;
  CScalingCalcOne(self->cscaling, x, NULL, &xc);

  int n_ls; DVRGetLSSize(self, &n_ls);
  PetscScalar *vs;  PetscMalloc1(n_ls, &vs);
  PetscInt    *idx; PetscMalloc1(n_ls, &idx);

  PetscBool zeroq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int i=0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int iq = 0; iq < nq; iq++) {
      idx[i] = i;
      ValueLS(self->bps, self->xs_c, nq, i_ele, iq, x, xc, &vs[i], &zeroq);
      i++;
    }
  }  

  Vec val_ls;
  DVRCreateR1LSVec(self, &val_ls);
  VecSetValues(val_ls, n_ls, idx, vs, INSERT_VALUES);
  VecAssemblyBegin(val_ls);
  VecAssemblyEnd(val_ls);

  Vec tmp;
  DVRCreateR1Vec(self, &tmp);
  MatMult(self->TT, val_ls, tmp);
  
  VecTDot(tmp, cs, y);
  
  PetscFree(vs);
  PetscFree(idx);
  VecDestroy(&tmp);
  VecDestroy(&val_ls);

  return 0;

}
PetscErrorCode DVRPsi(DVR self, Vec cs, Vec xs, Vec ys) {
  PetscErrorCode ierr;

  /* 
     yj = Sum(i)[ ci phi_i(xj) ] = Sum(i,k)[ ci Tik fk(xj)]
  */

  int nc; VecGetSize(cs, &nc);
  int nx; VecGetSize(xs, &nx);
  int ny; VecGetSize(ys, &ny);  
  if(nx != ny) {
    SETERRQ(self->comm, 1, "xs and ys must be same size");
  }
  if(self->num_basis != nc) {
    SETERRQ(self->comm, 1, "size of cs must be same as basis size");
  }

  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int nls; DVRGetLSSize(self, &nls);

  Mat fmat; MatCreate(self->comm, &fmat);
  MatSetSizes(fmat, PETSC_DECIDE, PETSC_DECIDE, nx, nls);
  MatSetUp(fmat);

  PetscScalar *ptr_x; VecGetArray(xs, &ptr_x);

  for(int ie = 0; ie < ne; ie++) {
    for(int iq = 0; iq < nq; iq++) {
      int kls = ie*nq + iq;
      for(int jx = 0; jx < nx; jx++) {
	PetscScalar y;
	PetscBool zeroq;
	ValueLS(self->bps, self->xs_c, self->nq,
		ie, iq, creal(ptr_x[jx]), ptr_x[jx], &y, &zeroq);
	if(!zeroq) {
	  MatSetValue(fmat, jx, kls, y, INSERT_VALUES);
	}
      }
    }
  }
  MatAssemblyBegin(fmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(fmat,   MAT_FINAL_ASSEMBLY);

  Vec TT_cs;
  ierr = MatCreateVecs(self->T, NULL, &TT_cs); CHKERRQ(ierr);
  ierr = MatMult(self->T, cs, TT_cs); CHKERRQ(ierr);
  ierr = MatMult(fmat, TT_cs, ys); CHKERRQ(ierr);

  MatDestroy(&fmat);
  VecDestroy(&TT_cs);
  VecRestoreArray(xs, &ptr_x);

  return 0;    
}
PetscErrorCode DVRDerivPsiOne(DVR self, Vec c, PetscReal x, PetscScalar *y) {

  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int nx = nq * ne;

  PetscScalar *vs;
  int *idx;
  PetscReal *zs; BPSGetZs(self->bps, &zs, NULL);
  PetscMalloc1(nx, &vs);
  PetscMalloc1(nx, &idx);

  PetscScalar xc; CScalingCalcOne(self->cscaling, x, NULL, &xc);

  int i=0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int iq = 0; iq < nq; iq++) {
      idx[i] = i;
      PetscBool is_in_q;
      ValueDerivLS(self->bps, self->xs_c, nq, i_ele, iq, x, xc, &vs[i], &is_in_q);
      i++;
    }
  }  

  Vec val_ls;
  DVRCreateR1LSVec(self, &val_ls);
  VecSetValues(val_ls, nx, idx, vs, INSERT_VALUES);
  VecAssemblyBegin(val_ls);
  VecAssemblyEnd(val_ls);

  Vec tmp;
  DVRCreateR1Vec(self, &tmp);
  MatMult(self->TT, val_ls, tmp);
  
  VecTDot(tmp, c, y);

  PetscFree(vs);
  PetscFree(idx);
  PetscFree(zs);
  VecDestroy(&val_ls);
  VecDestroy(&tmp);

  return 0;
}
PetscErrorCode DVRDerivPsi(DVR self, Vec cs, Vec xs, Vec ys) {
  PetscErrorCode ierr;

  /* 
     y'j = Sum(i)[ ci phi_i(xj) ] = Sum(i,k)[ ci Tik f'k(xj)]
  */

  int nc; VecGetSize(cs, &nc);
  int nx; VecGetSize(xs, &nx);
  int ny; VecGetSize(ys, &ny);  
  if(nx != ny) {
    SETERRQ(self->comm, 1, "xs and ys must be same size");
  }
  if(self->num_basis != nc) {
    SETERRQ(self->comm, 1, "size of cs must be same as basis size");
  }

  int ne; BPSGetNumEle(self->bps, &ne);
  int nq = self->nq;
  int nls; DVRGetLSSize(self, &nls);

  Mat fmat; MatCreate(self->comm, &fmat);
  MatSetSizes(fmat, PETSC_DECIDE, PETSC_DECIDE, nx, nls);
  MatSetUp(fmat);

  PetscScalar *ptr_x; VecGetArray(xs, &ptr_x);

  for(int ie = 0; ie < ne; ie++) {
    for(int iq = 0; iq < nq; iq++) {
      int kls = ie*nq + iq;
      for(int jx = 0; jx < nx; jx++) {
	PetscScalar y;
	PetscBool zeroq;
	ValueDerivLS(self->bps, self->xs_c, self->nq,
		     ie, iq, creal(ptr_x[jx]), ptr_x[jx], &y, &zeroq);
	if(!zeroq) {
	  MatSetValue(fmat, jx, kls, y, INSERT_VALUES);
	}
      }
    }
  }
  MatAssemblyBegin(fmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(fmat,   MAT_FINAL_ASSEMBLY);

  Vec TT_cs;
  ierr = MatCreateVecs(self->T, NULL, &TT_cs); CHKERRQ(ierr);
  ierr = MatMult(self->T, cs, TT_cs); CHKERRQ(ierr);
  ierr = MatMult(fmat, TT_cs, ys); CHKERRQ(ierr);

  MatDestroy(&fmat);
  VecDestroy(&TT_cs);
  VecRestoreArray(xs, &ptr_x);

  return 0;    
}

// ---- Inner function ----
PetscErrorCode DVRPrepareT2(DVR self) {
  PetscErrorCode ierr;
  ierr = MatMatSynthesize(self->T, self->T, 1.0, MAT_INITIAL_MATRIX, &self->T2); 
  CHKERRQ(ierr);  
  ierr = MatTranspose(self->T2, MAT_INITIAL_MATRIX, &self->T2T); 
  CHKERRQ(ierr);  
  return 0;
}

// ---- R1Mat/R2Mat ----
PetscErrorCode ScalarPower(MPI_Comm comm,
			   PetscScalar x, int n, PetscScalar* y) {
  if(n < 0) {
    SETERRQ(comm, 1, "n is negative");
  } else if(n == 0) {
    *y = 1.0;
  } else {

    *y = 1.0;
    for(int i = 0; i < n; i++) {
      *y *= x;
    }
  }

  return 0;
}
PetscErrorCode DVRCreateR1Vec(DVR self, Vec *m) {
  int n; DVRGetSize(self, &n);
  VecCreate(self->comm, m);
  VecSetSizes(*m, PETSC_DECIDE, n);
  VecSetUp(*m);
  return 0;
}
PetscErrorCode DVRCreateR1Mat(DVR self, Mat *M) {
  //  int ne; BPSGetNumEle(self->bps, &ne);
  //  int nq = self->nq;
  //  int n = (nq-2) * ne + ne -1;
  int n; DVRGetSize(self, &n);
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(*M);
  return 0;
}
PetscErrorCode DVRCreateR2Mat(DVR self, Mat *M) {
  int n; DVRGetSize(self, &n);
  MatCreate(self->comm, M);
  MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, n*n, n*n);
  MatSetUp(*M);
  return 0;
}

PetscErrorCode DVRPotR1Vec(DVR self, Pot pot, Vec v) {

  PetscErrorCode ierr;
  Vec LSR1;
  ierr = DVRCreateR1LSVec(self, &LSR1);  CHKERRQ(ierr);
  ierr = DVRPotR1LSVec(self, pot, LSR1); CHKERRQ(ierr);
  ierr = DVRLSVecToVec(self, LSR1, v);  CHKERRQ(ierr);
  ierr = VecDestroy(&LSR1);              CHKERRQ(ierr);
  return 0;  

  /*
  PetscErrorCode ierr;
  int n = self->num_basis;
  PetscScalar *vs;
  PetscMalloc1(n, &vs);
  ierr = PFApply(pot, n, self->xs_basis_c, vs);

  for(int i = 0; i < n; i++) {
    //    PetscScalar wi = self->ws_c[i];    
    PetscScalar wi = self->ws_basis_c[i];
    VecSetValue(v, i, csqrt(wi)*vs[i], INSERT_VALUES);
    //    printf("%f, %f, %f, %f\n", creal(wi), cimag(wi), creal(vs[i]), cimag(vs[i]));
  }

  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  PetscFree(vs);

  return 0;
  */
}
PetscErrorCode DVRSR1Mat(DVR self, Mat M) {
  PetscErrorCode ierr = 0;
  
  for(int i = 0; i < self->num_basis; i++) 
    ierr = MatSetValue(M, i, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode DVRD2R1Mat(DVR self, Mat M) {
  PetscErrorCode ierr;
  Mat LSR1, T;
  ierr = DVRCreateR1LSMat(self, &LSR1); CHKERRQ(ierr);
  //  DVRCreateR1Mat(self, &T);
  ierr = DVRD2R1LSMat(self, LSR1); CHKERRQ(ierr);
  ierr = DVRLSMatToMat(self, LSR1, MAT_INITIAL_MATRIX, &T); CHKERRQ(ierr);
  ierr = MatDestroy(&LSR1); CHKERRQ(ierr);
  ierr = MatCopy(T, M, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatDestroy(&T); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode DVRPotR1Mat(DVR self, Pot pot, Mat M) {
  PetscErrorCode ierr;

  PetscScalar *vs;
  
  PetscMalloc1(self->num_basis, &vs);
  ierr = PFApply(pot, self->num_basis, self->xs_basis_c, vs);
  
  for(int i = 0; i < self->num_basis; i++) {
    ierr = MatSetValue(M, i, i, vs[i], INSERT_VALUES); CHKERRQ(ierr);
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  PetscFree(vs);
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
    PetscScalar xc = self->xs_basis_c[i];
    PetscScalar g = x > a ? xc : a;
    PetscScalar s = x < a ? xc : a;
    PetscScalar v=0.0;
    //ScalarPower(self->comm, q, s/g, &v);
    ScalarPower(self->comm, s/g, q, &v);
    v /= g;
    ierr = MatSetValue(M, i, i, v, INSERT_VALUES); CHKERRQ(ierr);

  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVREER2Mat_direct(DVR self, int q, Mat M) {
 PetscErrorCode ierr;
  int num; ierr = DVRGetSize(self, &num); CHKERRQ(ierr);
  for(int i = 0; i < num; i++) {
    for(int j =  0; j < num; j++) {
      PetscScalar ri = self->xs_basis_c[i];
      PetscScalar rj = self->xs_basis_c[j];
      PetscScalar g,s;
      if(PetscRealPart(ri) > PetscRealPart(rj)) {
	g = ri; s = rj;
      } else {
	g = rj; s = ri;
      }
      PetscScalar val=0.0;
      ScalarPower(self->comm, q, s/g, &val);
      val /= g;
      int idx = i*num+j;
      ierr = MatSetValue(M, idx, idx, val, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  M, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode DVREER2Mat_poisson(DVR self, int q, Mat M) {
  PetscErrorCode ierr;
  Mat T;
  ierr = DVRCreateR1Mat(self, &T); CHKERRQ(ierr);
  ierr = DVRD2R1Mat(self, T); CHKERRQ(ierr);
  ierr = MatScale(T, -1.0); CHKERRQ(ierr);
  if(q != 0) {
    Mat L2;
    ierr = DVRCreateR1Mat(self, &L2); CHKERRQ(ierr);
    ierr = DVRR2invR1Mat(self, L2); CHKERRQ(ierr);
    ierr = MatAXPY(T, q*(q+1.0), L2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    MatDestroy(&L2);
  }

  int num; ierr = DVRGetSize(self, &num); CHKERRQ(ierr);
  Vec c; ierr = DVRCreateR1Vec(self, &c); CHKERRQ(ierr);
  Vec y; ierr = DVRCreateR1Vec(self, &y); CHKERRQ(ierr);
  PetscScalar *ys; ierr = PetscMalloc1(num, &ys); CHKERRQ(ierr);
  KSP ksp;
  ierr = KSPCreate(self->comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, T, T); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  for(int i = 0; i < num; i++) {
    // -- Set unit vector --
    for(int j = 0; j < num; j++) {
      ierr = VecSetValue(c, j, i==j?1.0:0.0, INSERT_VALUES); CHKERRQ(ierr);
    }
    VecAssemblyBegin(c); VecAssemblyEnd(c);

    // -- Solve linear problem --
    ierr = KSPSolve(ksp, c, y); CHKERRQ(ierr);
    ierr = VecGetArray(y, &ys); CHKERRQ(ierr);

    // -- compute ERI --
    for(int j =  0; j < num; j++) {
      int idx = i*num+j;
      PetscScalar ri = self->xs_basis_c[i];
      PetscScalar rj = self->xs_basis_c[j];
      PetscScalar val = ys[j]*(2*q+1) * self->vs_basis[j] * self->vs_basis[i] / (ri*rj);
      PetscScalar rmax = self->xs_basis_c[num-1];
      val += pow(ri * rj, q) / pow(rmax, 2*q+1);
      ierr = MatSetValue(M, idx, idx, val, INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(y, &ys); CHKERRQ(ierr);
   
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

  // -- Finalize --
  MatDestroy(&T);  
  PetscFree(ys);
  VecDestroy(&c);
  VecDestroy(&y);
  KSPDestroy(&ksp);
  return 0;
}
PetscErrorCode DVREER2Mat(DVR self, int q, Mat M) {

  PetscErrorCode ierr;
  //  DVREER2Mat_direct(self, q, M);
  ierr = DVREER2Mat_poisson(self, q, M); CHKERRQ(ierr);
  return 0;

  /*
  PetscErrorCode ierr;
  Mat LS;
  ierr = DVRCreateR1LSMat(self, &LS); CHKERRQ(ierr);
  ierr = DVREER2LSMat(self, q, LS); CHKERRQ(ierr);
  ierr = DVRR2LSMatToR2Mat(self, LS, &M); CHKERRQ(ierr);
  return 0;
  */
}

// ---- LSR1Mat/LSR2Mat ----
PetscErrorCode DVRCreateR1LSVec(DVR self, Vec *m) {

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int n = nq * ne;
  VecCreate(self->comm, m);
  VecSetSizes(*m, PETSC_DECIDE, n);
  VecSetUp(*m);
  return 0;

}
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

PetscErrorCode DVRPotR1LSVec(DVR self, Pot pot, Vec v) {

  PetscErrorCode ierr;

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int n = nq * ne;
  PetscScalar *vs; PetscMalloc1(n, &vs);    
  ierr = PFApply(pot, n, self->xs_c, vs);

  for(int k = 0; k < n; k++)
    VecSetValue(v, k, self->ws_c[k]*vs[k], INSERT_VALUES);

  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  PetscFree(vs);
  
  return 0;
}
PetscErrorCode DVRSR1LSMat(DVR self, Mat M) {

  int nq = self->nq;
  int ne; BPSGetNumEle(self->bps, &ne);
  int n = nq * ne;

  for(int k = 0; k < n; k++)
    MatSetValue(M, k, k, self->ws_c[k], INSERT_VALUES);

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
	DerivLS(self->xs_c, self->ws_c,
		ne, nq, i, m, mp, &ds_list[idx]);
      }

  DVRCreateR1LSMat(self, &M);

  for(int i = 0; i < ne; i++) 
    for(int ma = 0; ma < nq; ma++)
      for(int mb = 0; mb < nq; mb++) {
	int idx_a = i*nq + ma;
	int idx_b = i*nq + mb;
	PetscScalar v = 0.0;
	for(int k = 0; k < nq; k++)
	  v -= ds_list[idx_a*nq+k] * ds_list[idx_b*nq+k] * self->ws_c[i*nq+k];
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
      x = self->xs_c[idx];
      w = self->ws_c[idx];
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

  int numi, numj;
  MatGetSize(M, &numi, &numj);
  if(numi != nq*ne || numj != nq*ne) {
    SETERRQ(self->comm, 1, "Size mismatch.");
  }

  int idx = 0;
  for(int i_ele = 0; i_ele < ne; i_ele++) {
    for(int i_q = 0; i_q < nq;   i_q++) {
      PetscReal x = x = self->xs[idx];
      PetscScalar xc = self->xs_c[idx];
      PetscScalar wc = self->ws_c[idx];
      PetscScalar v = 0.0;
      if(i_ele == 0 && i_q == 0)
	v = 0.0;
      else if(q == 0) {
	v = 1.0/xc;
      }
      else {
	PetscScalar g = x > a ? xc : a;
	PetscScalar s = x < a ? xc : a;
	ScalarPower(self->comm, s/g, q, &v);
	v /= g;
      }
      ierr = MatSetValue(M, idx, idx, v*wc, INSERT_VALUES);
      CHKERRQ(ierr);
      idx++;
    }    
  }  

  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  
  return 0;
}

PetscErrorCode DVRLSVecToVec(DVR self, Vec A, Vec B) {

  PetscErrorCode ierr;
  ierr = MatMult(self->TT, A, B); CHKERRQ(ierr);
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
	  double ri = self->xs[i_ele*nq+i]; double wi = self->ws_c[i_ele*nq+i];
	  double rj = self->xs[j_ele*nq+j]; double wj = self->ws_c[j_ele*nq+j];
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

