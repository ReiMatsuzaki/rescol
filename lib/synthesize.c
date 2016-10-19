#include "../include/synthesize.h"

PetscErrorCode VecVecSynthesizeSymbolic(Vec A, Vec B, Vec *C) {
  /*
    perform construction, preallocation, and computes the ij structure
    see MatMatMultSymbolc for example.
   */

  PetscErrorCode ierr;
  MPI_Comm comm;   ierr = PetscObjectGetComm((PetscObject)A, &comm);

  // construction
  PetscInt na; VecGetSize(A, &na);
  PetscInt nb; VecGetSize(B, &nb);
  ierr =VecCreate(comm, C);  CHKERRQ(ierr);
  ierr =VecSetSizes(*C, PETSC_DECIDE, na*nb);  CHKERRQ(ierr);
  ierr =VecSetUp(*C);  CHKERRQ(ierr);

  // ij structure (values are all zeo)
  PetscScalar *as; VecGetArray(A, &as);
  PetscScalar *bs; VecGetArray(B, &bs);
  PetscScalar *cs; PetscMalloc1(na*nb, &cs);
  PetscInt *idxs; PetscMalloc1(na*nb, &idxs);
  
  PetscInt idx = 0;
  for(int j = 0; j < nb; j++) 
    for(int i = 0; i < na; i++) {
      cs[i + na*j] = 0.0;
      idxs[i + na*j] = idx;
      idx++;
    }
  ierr = VecSetValues(*C, na*nb, idxs, cs, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(*C);  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*C); CHKERRQ(ierr);

  VecRestoreArray(A, &as); 
  VecRestoreArray(B, &bs); 
  PetscFree(cs);
  PetscFree(idxs);  
  
  return 0;
}
PetscErrorCode VecVecSynthesizeNumeric(Vec A, Vec B, PetscScalar c, Vec C) {

  PetscInt na; VecGetSize(A, &na);
  PetscInt nb; VecGetSize(B, &nb);

  PetscScalar *as; VecGetArray(A, &as);
  PetscScalar *bs; VecGetArray(B, &bs);
  PetscScalar *cs; VecGetArray(C, &cs);

  for(int j = 0; j < nb; j++) {
    PetscScalar bc = bs[j] * c;
    for(int i = 0; i < na; i++) {
      cs[i + na*j] = as[i] * bc;
    }
  }

  VecRestoreArray(A, &as); 
  VecRestoreArray(B, &bs); 
  VecRestoreArray(C, &cs); 
  return 0;
}
PetscErrorCode VecVecSynthesize(Vec A, Vec B, PetscScalar c, 
			     MatReuse scall, Vec *C){
  PetscErrorCode ierr;
  MPI_Comm comm;   ierr = PetscObjectGetComm((PetscObject)A, &comm);

  if(scall == MAT_INITIAL_MATRIX) {
    ierr = VecVecSynthesizeSymbolic(A, B, C); CHKERRQ(ierr);
  } else if(scall == MAT_REUSE_MATRIX) {

  } else
    SETERRQ(comm, 1, "MatReuse scall <- {MAT_INITIAL_MATRIX, MAT_REUSE_MATRIX}");
    
  ierr = VecVecSynthesizeNumeric(A, B, c, *C); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MatMatSynthesizeSymbolic_Default(Mat A, Mat B, Mat *C) {
  PetscErrorCode ierr;
  MPI_Comm comm;   PetscObjectGetComm((PetscObject)A, &comm);

  // construction  
  PetscInt na, ma; ierr = MatGetSize(A, &na, &ma); CHKERRQ(ierr);
  PetscInt nb, mb; ierr = MatGetSize(B, &nb, &mb); CHKERRQ(ierr);
  
  ierr = MatCreate(comm, C); CHKERRQ(ierr);
  ierr = MatSetSizes(*C, PETSC_DECIDE, PETSC_DECIDE, na*nb, ma*mb); CHKERRQ(ierr);
  ierr = MatSetUp(*C); CHKERRQ(ierr);  

  return 0;
}
PetscErrorCode MatMatSynthesizeNumeric_Default(Mat A, Mat B, PetscScalar c, Mat C) {
  PetscErrorCode ierr;
  MPI_Comm comm;   PetscObjectGetComm((PetscObject)A, &comm);
  PetscInt na, ma; ierr = MatGetSize(A, &na, &ma); CHKERRQ(ierr);
  PetscInt nb, mb; ierr = MatGetSize(B, &nb, &mb); CHKERRQ(ierr);
  PetscInt nc, mc; ierr = MatGetSize(C, &nc, &mc); CHKERRQ(ierr);
  
  const PetscScalar **row_a, **row_b;
  PetscInt *ncols_a, *ncols_b;
  const PetscInt **cols_a, **cols_b;
  
  PetscMalloc1(na, &row_a);   PetscMalloc1(nb, &row_b);
  PetscMalloc1(na, &ncols_a); PetscMalloc1(nb, &ncols_b);
  PetscMalloc1(na, &cols_a);  PetscMalloc1(nb, &cols_b);

  for(int i = 0; i < na; i++) 
    ierr = MatGetRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  for(int i = 0; i < nb; i++) 
    ierr = MatGetRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);

  for(int i_a = 0; i_a < na; i_a++) {
    for(int idx_a = 0; idx_a < ncols_a[i_a]; idx_a++) {
      int j_a = cols_a[i_a][idx_a];
      for(int i_b = 0; i_b < nb; i_b++) { 
      	for(int idx_b = 0; idx_b < ncols_b[i_b]; idx_b++) {
	  int j_b = cols_b[i_b][idx_b];
	  int i = i_a + i_b * na;
	  int j = j_a + j_b * ma;
	  PetscScalar v =  row_a[i_a][idx_a] * row_b[i_b][idx_b] * c;
	  ierr = MatSetValue(C, i, j, v, INSERT_VALUES);  CHKERRQ(ierr);
	}
      }
    }
  }

  for(int i = 0; i < na; i++)
    ierr = MatRestoreRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  for(int i = 0; i < nb; i++)
    ierr = MatRestoreRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
  PetscFree(row_a);   PetscFree(row_b);
  PetscFree(ncols_a); PetscFree(ncols_b);
  PetscFree(cols_a);  PetscFree(cols_b);

  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode MatMatSynthesizeSymbolic_Fast(Mat A, Mat B, Mat *C) {

  MPI_Comm comm;   PetscObjectGetComm((PetscObject)A, &comm);

  PetscInt an, am, bn, bm;
  MatGetSize(A, &am, &an);
  MatGetSize(B, &bm, &bn);

  PetscInt cm = am*bm;
  PetscInt cn = an*bn;
  
  MatCreate(comm, C);
  MatSetSizes(*C, PETSC_DECIDE, PETSC_DECIDE, cm, cn);
  MatSetUp(*C);
  
  return 0;
}
PetscErrorCode MatMatSynthesizeNumeric_Fast(Mat A, Mat B, PetscScalar c, Mat C) {
  PetscErrorCode ierr;

  PetscInt an, am, bn, bm;
  MatGetSize(A, &am, &an);
  MatGetSize(B, &bm, &bn);

  const PetscScalar **row_a, **row_b;
  PetscInt *ncols_a, *ncols_b;
  const PetscInt **cols_a, **cols_b;
  
  PetscMalloc1(am, &row_a);   PetscMalloc1(bm, &row_b);
  PetscMalloc1(am, &ncols_a); PetscMalloc1(bm, &ncols_b);
  PetscMalloc1(am, &cols_a);  PetscMalloc1(bm, &cols_b);

  for(int i = 0; i < am; i++) {
    ierr = MatGetRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  }
  for(int i = 0; i < bm; i++) {
    ierr = MatGetRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
  }

  for(int ia = 0; ia < am; ia++) 
    for(int ib = 0; ib < bm; ib++) {
      PetscInt ncols = ncols_a[ia] * ncols_b[ib];
      PetscScalar *ac; PetscMalloc1(ncols, &ac);      
      PetscInt *idxn; PetscMalloc1(ncols, &idxn);
      int idxm[1] = {ia + am*ib};
      for(int idx_a = 0; idx_a < ncols_a[ia]; idx_a++)
	for(int idx_b = 0; idx_b < ncols_b[ib]; idx_b++) {
	  int ja = cols_a[ia][idx_a];
	  int jb = cols_b[ib][idx_b];
	  int idx = idx_a + idx_b*ncols_a[ia];
	  idxn[idx] = ja + jb * an;
	  ac[idx] = row_a[ia][idx_a] * row_b[ib][idx_b] * c;
	}
      MatSetValues(C, 1, idxm, ncols, idxn, ac, INSERT_VALUES);
      PetscFree(ac);
      PetscFree(idxn);
    }

  for(int i = 0; i < am; i++)
    ierr = MatRestoreRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  for(int i = 0; i < bm; i++)
    ierr = MatRestoreRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
  PetscFree(row_a);   PetscFree(row_b);
  PetscFree(ncols_a); PetscFree(ncols_b);
  PetscFree(cols_a);  PetscFree(cols_b);

  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);  

  return 0;
}
PetscErrorCode MatMatSynthesizeNumeric_Fast2(Mat A, Mat B, PetscScalar c, Mat C) {
  PetscErrorCode ierr;

  PetscInt an, am, bn, bm;
  MatGetSize(A, &am, &an);
  MatGetSize(B, &bm, &bn);

  const PetscScalar **row_a, **row_b;
  PetscInt *ncols_a, *ncols_b;
  const PetscInt **cols_a, **cols_b;
  
  PetscMalloc1(am, &row_a);   PetscMalloc1(bm, &row_b);
  PetscMalloc1(am, &ncols_a); PetscMalloc1(bm, &ncols_b);
  PetscMalloc1(am, &cols_a);  PetscMalloc1(bm, &cols_b);

  int ncol_a_max = 0, ncol_b_max = 0;
  for(int i = 0; i < am; i++) {
    ierr = MatGetRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
    if(ncol_a_max < ncols_a[i])
      ncol_a_max = ncols_a[i];
  }
  for(int i = 0; i < bm; i++) {
    ierr = MatGetRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
    if(ncol_b_max < ncols_b[i])
      ncol_b_max = ncols_b[i];
  }

  PetscScalar *ac; PetscMalloc1(ncol_a_max*ncol_b_max, &ac);
  PetscInt *idxn; PetscMalloc1(ncol_a_max*ncol_b_max, &idxn);
  PetscScalar *space; PetscMalloc1(ncol_b_max*bm, &space);

  for(int ib = 0; ib < bm; ib++) {
    for(int idx_b = 0; idx_b < ncols_b[ib]; idx_b++) {
      space[ib+bm*idx_b] = row_b[ib][idx_b] * c;
    }
  }

  for(int ia = 0; ia < am; ia++) 
    for(int ib = 0; ib < bm; ib++) {
      int idxm[1] = {ia + am*ib};
      PetscInt ncols = ncols_a[ia] * ncols_b[ib];

      for(int idx_a = 0; idx_a < ncols_a[ia]; idx_a++) {
	PetscScalar r = row_a[ia][idx_a];
	for(int idx_b = 0; idx_b < ncols_b[ib]; idx_b++) {
	  int ja = cols_a[ia][idx_a];
	  int jb = cols_b[ib][idx_b];
	  int idx = idx_a + idx_b*ncols_a[ia];
	  idxn[idx] = ja + jb * an;
	  ac[idx] = r * space[ib+bm*idx_b];
	}
      }
      MatSetValues(C, 1, idxm, ncols, idxn, ac, INSERT_VALUES);
    }

  for(int i = 0; i < am; i++)
    ierr = MatRestoreRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  for(int i = 0; i < bm; i++)
    ierr = MatRestoreRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
  PetscFree(row_a);   PetscFree(row_b);
  PetscFree(ncols_a); PetscFree(ncols_b);
  PetscFree(cols_a);  PetscFree(cols_b);
  
  PetscFree(ac); PetscFree(idxn); PetscFree(space);

  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);  

  return 0;
}
PetscErrorCode MatMatSynthesizeSymbolic_SeqAIJ(Mat A, Mat B, Mat *C) {
  SETERRQ(PetscObjectComm((PetscObject)A), 1, "not implemented yet");

/*
  PetscErrorCode     ierr;
  Mat_SeqAIJ         *a =(Mat_SeqAIJ*)A->data,*b=(Mat_SeqAIJ*)B->data,*c;
  PetscInt           *ai=a->i,*bi=b->i, *aj=a->j, *bj=b->j, *ci,*cj;
  PetscScalar *ca;

  // A is (am, an) matrix; B is(bm, bn) matrix
  PetscInt           am =A->rmap->N,bm=B->rmap->N;
  PetscInt           an =A->cmap->N;//,bn=B->cmap->N;

  ierr = PetscMalloc1(am*an, &ci); CHKERRQ(ierr);
  ierr = PetscMalloc1(ai[am]*bi[bm], &cj); CHKERRQ(ierr);
  ierr = PetscMalloc1(ai[am]*bi[bm]+1, &ca); CHKERRQ(ierr);
  ierr = PetscMemzero(ca, (ai[am]*bi[bm]+1)*sizeof(PetscScalar)); CHKERRQ(ierr);

  int idx=0;
  for(int k_a = 0; k_a < am; k_a++)
    for(int k_b = 0; k_b < bm; k_b++) {
      int k_c = k_a + k_b*am;
      ci[k_c] = idx;
      for(int idx_a = ai[k_a]; idx_a < ai[k_a+1]; idx_a++)
	for(int idx_b = bi[k_b]; idx_b < ai[k_b+1]; idx_b++) {
	  cj[idx] = aj[idx_a-1] + bj[idx_b-1]*an;
	  idx++;
	}
    }
  ci[am*bm-1] = idx;

  ierr = MatCreateSeqAIJWithArrays(PetscObjectComm((PetscObject)A), 
				   am*bm,an*am,ci,cj,ca,C);
  c  = (Mat_SeqAIJ*)((*C)->data);
  c->free_a  = PETSC_TRUE;
  c->free_ij = PETSC_TRUE;
  c->nonew   = 0;
*/
  return 0;
}
PetscErrorCode MatMatSynthesizeNumeric_SeqAIJ(Mat A, Mat B, PetscScalar x, Mat C) {
  SETERRQ(PetscObjectComm((PetscObject)A), 1, "not implemented yet");
/*
  PetscErrorCode ierr;
  Mat_SeqAIJ     *a =(Mat_SeqAIJ*)A->data,*b=(Mat_SeqAIJ*)B->data, *c;
  PetscInt       *ai=a->i,*bi=b->i, *aj=a->j, *bj=b->j, *ci,*cj;
  //  PetscScalar    *aa, *ba, *ca;

  // A is (am, an) matrix; B is(bm, bn) matrix
  PetscInt           am =A->rmap->N,bm=B->rmap->N;
  PetscInt           an =A->cmap->N;//,bn=B->cmap->N;

  ierr = PetscMalloc1(am*an, &ci); CHKERRQ(ierr);
  ierr = PetscMalloc1(ai[am]*bi[bm], &cj); CHKERRQ(ierr);
  ierr = PetscMalloc1(ai[am]*bi[bm]+1, &ca); CHKERRQ(ierr);
  ierr = PetscMemzero(ca, (ai[am]*bi[bm]+1)*sizeof(PetscScalar)); CHKERRQ(ierr);

  int idx=0;
  for(int k_a = 0; k_a < am; k_a++)
    for(int k_b = 0; k_b < bm; k_b++) {
      int k_c = k_a + k_b*am;
      ci[k_c] = idx;
      for(int idx_a = ai[k_a]; idx_a < ai[k_a+1]; idx_a++)
	for(int idx_b = bi[k_b]; idx_b < ai[k_b+1]; idx_b++) {
	  cj[idx] = aj[idx_a-1] + bj[idx_b-1]*an;
	  c[idx] = a[idx_a-1] * b[idx_b-1] * c;
	  idx++;
	}
    }
  ci[am*bm-1] = idx;

  ierr = MatCreateSeqAIJWithArrays(PetscObjectComm((PetscObject)A), 
				   am*bm,an*am,ci,cj,ca,C);
  c  = (Mat_SeqAIJ*)((*C)->data);
  c->free_a  = PETSC_TRUE;
  c->free_ij = PETSC_TRUE;
  c->nonew   = 0;
*/
  return 0;
}
PetscErrorCode MatMatSynthesizeSymbolic(Mat A, Mat B, Mat *C) {
  PetscErrorCode ierr;
  int alg = 0;
  PetscOptionsGetInt(NULL, NULL, "-mat_mat_synthesize_alg", &alg, NULL);
  switch(alg) {
  case 1:
    ierr = MatMatSynthesizeSymbolic_SeqAIJ(A, B, C); CHKERRQ(ierr);
    break;
  case 2:
    ierr = MatMatSynthesizeSymbolic_Fast(A, B, C); CHKERRQ(ierr);
    break;
  case 3:
    ierr = MatMatSynthesizeSymbolic_Fast(A, B, C); CHKERRQ(ierr);
      break;
  default:
    ierr = MatMatSynthesizeSymbolic_Default(A, B, C); CHKERRQ(ierr);
    break;
  }

  return 0;
}
PetscErrorCode MatMatSynthesizeNumeric(Mat A, Mat B, PetscScalar a, Mat C) {
  PetscErrorCode ierr;
  int alg = 0;
  PetscOptionsGetInt(NULL, NULL, "-mat_mat_synthesize_alg", &alg, NULL);

  switch(alg) {
  case 1:
    ierr = MatMatSynthesizeNumeric_SeqAIJ(A, B, a, C); CHKERRQ(ierr);    
    break;
  case 2:
    ierr = MatMatSynthesizeNumeric_Fast(A, B, a, C); CHKERRQ(ierr);    
    break;
  case 3:
    ierr = MatMatSynthesizeNumeric_Fast2(A, B, a, C); CHKERRQ(ierr);    
    break;
  default:
    ierr = MatMatSynthesizeNumeric_Default(A, B, a, C); CHKERRQ(ierr);    
    break;
  }
  return 0;
}
PetscErrorCode MatMatSynthesize(Mat A, Mat B, PetscScalar a, MatReuse scall, Mat *C){
  PetscErrorCode ierr;
  if(scall == MAT_INITIAL_MATRIX) {
    ierr = MatMatSynthesizeSymbolic(A, B, C); CHKERRQ(ierr);
  }
  ierr = MatMatSynthesizeNumeric(A, B, a, *C); CHKERRQ(ierr);
 return 0;

}

PetscErrorCode MatMatMatSynthesizeSymbolic(Mat A, Mat B, Mat C, Mat *D) {

  PetscErrorCode ierr;
  Mat BC; 
  ierr = MatMatSynthesizeSymbolic(B, C, &BC); CHKERRQ(ierr);
  ierr = MatMatSynthesizeSymbolic(A, BC, D); CHKERRQ(ierr);
  MatDestroy(&BC);
  return 0;

}
PetscErrorCode MatMatMatSynthesizeNumeric(Mat A, Mat B, Mat C, PetscScalar d, Mat D) {

  PetscErrorCode ierr;
  Mat BC; 
  ierr = MatMatSynthesizeSymbolic(B, C, &BC); CHKERRQ(ierr);
  ierr = MatMatSynthesizeNumeric(B, C, d, BC); CHKERRQ(ierr);
  ierr = MatMatSynthesizeNumeric(A, BC, 1.0, D); CHKERRQ(ierr);
  MatDestroy(&BC);
  return 0;
  
}
PetscErrorCode MatMatMatSynthesize(Mat A, Mat B, Mat C, 
				   PetscScalar d, MatReuse scall, Mat *D) {

  PetscErrorCode ierr;
  MPI_Comm comm;   ierr = PetscObjectGetComm((PetscObject)A, &comm);

  if(scall == MAT_INITIAL_MATRIX) {
    ierr = MatMatMatSynthesizeSymbolic(A, B, C, D); CHKERRQ(ierr);
  } else if(scall == MAT_REUSE_MATRIX) {

  } else
    SETERRQ(comm, 1, "MatReuse scall <- {MAT_INITIAL_MATRIX, MAT_REUSE_MATRIX}");
  ierr = MatMatMatSynthesizeNumeric(A, B, C, d, *D); CHKERRQ(ierr);
  return 0;
}

