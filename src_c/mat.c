#include "mat.h"

PetscErrorCode MatCreateFromCOOFormatFileOld(char* path, Mat* mat) {
  FILE* fp;
  PetscInt *rows, *cols;
  PetscScalar *datas;
  PetscInt i;
  PetscErrorCode ierr;
  int num_data, num_row, num_col, ret;

  if((fp = fopen(path, "r")) == NULL) {
    const char* msg = "Failed to open file.\0";
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, msg);
  }
  
  ret = fscanf(fp, "%d %d %d", &num_data, &num_row, &num_col);
  if(ret == EOF) {
    const char* msg = "Failed to read first line. Expected format is:\n num_data, num_row, num_col\0";
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, msg);
  }
    
  rows = (PetscInt*)malloc(sizeof(PetscInt)*num_data);
  cols = (PetscInt*)malloc(sizeof(PetscInt)*num_data);
  datas = (PetscScalar*)malloc(sizeof(PetscScalar)*num_data);
  if(rows == NULL || cols == NULL || datas == NULL){
    const char* msg = "Failed to allocate memory for row or col or data\0";
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_MEM, msg);
  }

  i = 0;
  while((ret = fscanf(fp, "%d %d %lf", &rows[i], &cols[i], &datas[i])) != EOF) { i++; }

  ierr = MatCreateSeqAIJFromTriple(PETSC_COMM_WORLD, 
				   num_row, num_col, 
				   rows, cols, datas, 
				   mat,
				   num_data, 0); CHKERRQ(ierr);
  fclose(fp);
  return ierr;
}

PetscErrorCode MatCreateFromCOOFormatFileHandler(FILE* fp, Mat* mat) {

  PetscInt col, row;
  PetscErrorCode ierr;
  int num_data, num_row, num_col;
  PetscScalar dat;  

  if(fp == NULL) {
    char msg[256] = "file path is NULL";
    SETERRQ(PETSC_COMM_WORLD, 1, msg);
  }

  if(fscanf(fp, "%d %d %d", &num_row, &num_col, &num_data) == EOF) {
    const char* msg = "Failed to read first line. Expected format is:\n num_row, num_col, num_data\0";
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, msg);
  }
  
  if(num_data <= 0 || num_row <= 0 || num_col <= 0 || num_row != num_col) {
    char msg[256];
    sprintf(msg, "Invalid head value. num_data=%d, num_row=%d, num_col=%d", 
	    num_data, num_row, num_col);
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, msg);    
  }

  // Create Mat
  ierr = MatCreate(PETSC_COMM_WORLD, mat); CHKERRQ(ierr);
  ierr = MatSetSizes(*mat, PETSC_DECIDE, PETSC_DECIDE, num_row, num_col); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*mat); CHKERRQ(ierr);
  ierr = MatSetUp(*mat); CHKERRQ(ierr);

  while(fscanf(fp, "%d %d %lf", &row, &col, &dat) != EOF) {
    ierr = MatSetValue(*mat, row, col, dat, INSERT_VALUES); CHKERRQ(ierr);
  }
  MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);
  return 0;
}

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat) {

  PetscErrorCode ierr;
  FILE* fp = NULL;

  if((fp = fopen(path, "r")) == NULL) {
    char msg[256];
    sprintf(msg, "Failed to open file. Target file path is :%s", path);
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, msg);
  }

  ierr = MatCreateFromCOOFormatFileHandler(fp, mat); CHKERRQ(ierr);

  fclose(fp);
  return 0;
}

PetscErrorCode VecCreateFromFile(const char* path, MPI_Comm comm, Vec *v ) {

  PetscErrorCode ierr;
  FILE *fp = NULL;
  
  if((fp = fopen(path, "r")) == NULL) {
    char msg[256]; 
    sprintf(msg, "Failed to open file: %s", path);
    SETERRQ(comm, 1, "failed to open file");
  }

  PetscInt n;
  if(fscanf(fp, "%d", &n) == EOF) {
    SETERRQ(comm, 1, "failed to read header (# of data)");
  }

  VecCreate(comm, v);
  VecSetSizes(*v, PETSC_DECIDE, n);
  VecSetFromOptions(*v);
  VecSetUp(*v);

  PetscInt i = 0;
  double x;
  while(fscanf(fp, "%lf", &x) != EOF) {
    ierr = VecSetValue(*v, i, x, INSERT_VALUES); CHKERRQ(ierr);
    i++;
  }
  fclose(fp);
  VecAssemblyBegin(*v);VecAssemblyEnd(*v);

  return 0;
}

PetscErrorCode MatSetDirFile(const char* dn, const char* fn, Mat *M) {
  PetscErrorCode ierr;
  char path[100];
  sprintf(path, "%s/%s", dn, fn);
  ierr = MatCreateFromCOOFormatFile(path, M); CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode PrintTimeStamp(MPI_Comm comm, const char* label, time_t *t) {
  time_t tt;
  time(&tt);
  if(t != NULL)
    *t = tt;
  PetscPrintf(comm, "[%10s] %s", label, ctime(&tt));
  return 0;
}

PetscErrorCode EPSWriteToFile(EPS eps, char* path_detail, char* path_eigvals, char* path_eigvecs) {

  /*
    Parameters
    ----------
    eps : EPS Context
    path_detail : file path for writing calculation detail
    path_eigvals : file path for eigenvalus (if NULL, no output)
    path_eigvecs : file path for eigenvectors (if NULL, no output)
   */
  
  FILE* fp_detail = NULL;
  FILE* fp_eigvals = NULL;
  FILE* fp_eigvecs = NULL;
  EPSType type;
  PetscErrorCode ierr;
  PetscInt nconv, i, its, nev, maxit;
  PetscScalar eig, im_eig, error, tol;
  Mat H;
  Vec xs, ys;
  PetscViewer vec_viewer;

  // prepare vector
  ierr = EPSGetOperators(eps, &H, NULL); CHKERRQ(ierr);
  ierr = MatCreateVecs(H, NULL, &xs); CHKERRQ(ierr);
  ierr = MatCreateVecs(H, NULL, &ys); CHKERRQ(ierr);
  if(path_eigvecs)
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, path_eigvecs, &vec_viewer);

  // open file
  fp_detail = fopen(path_detail, "w");
  if(path_eigvals)
    fp_eigvals = fopen(path_eigvals, "w");
  if(path_eigvecs)
    fp_eigvecs = fopen(path_eigvecs, "w");

  // extract basic information
  EPSGetIterationNumber(eps, &its); 
  PetscFPrintf(PETSC_COMM_WORLD, fp_detail, "iteration: %d\n", its);
  EPSGetType(eps, &type); 
  PetscFPrintf(PETSC_COMM_WORLD, fp_detail, "EPSType: %s\n", type);
  EPSGetDimensions(eps, &nev, NULL, NULL);
  PetscFPrintf(PETSC_COMM_WORLD, fp_detail, "dim: %d\n", nev);
  ierr = EPSGetTolerances(eps, &tol, &maxit);
  PetscFPrintf(PETSC_COMM_WORLD, fp_detail, "tol: %f\n", tol);
  PetscFPrintf(PETSC_COMM_WORLD, fp_detail, "maxit: %d\n", maxit);

  // write to files
  PetscFPrintf(PETSC_COMM_WORLD, fp_eigvals, "EigenValue, Error\n");
  ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
  for(i = 0; i < nconv; i++) {
    ierr = EPSGetEigenpair(eps, i, &eig, &im_eig, xs, ys); CHKERRQ(ierr);
    ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error); CHKERRQ(ierr);
    if(fp_eigvals)
      PetscFPrintf(PETSC_COMM_WORLD, fp_eigvals, "%g %g\n", eig, error);    
    if(fp_eigvecs)
      VecView(xs, vec_viewer);
  }
  
  // finalize
  VecDestroy(&xs);
  VecDestroy(&ys);
  fclose(fp_detail);
  if(fp_eigvals)
    fclose(fp_eigvals);
  if(fp_eigvecs)
    fclose(fp_eigvecs);

  return 0;
}

PetscErrorCode VecInitSynthesize(Vec A, Vec B, MPI_Comm comm, Vec *C) {
  
  PetscInt na, nb;
  VecGetSize(A, &na); VecGetSize(B, &nb);

  VecCreate(comm, C);
  VecSetSizes(*C, PETSC_DECIDE, na*nb);
  VecSetFromOptions(*C);
  
  return 0;
}

PetscErrorCode VecSynthesize(Vec A, Vec B, PetscScalar c, 
			     Vec *C, InsertMode mode) {
  
  PetscInt na, nb;
  VecGetSize(A, &na); VecGetSize(B, &nb);

  PetscScalar *as, *bs, *cs;
  VecGetArray(A, &as); VecGetArray(B, &bs);
  cs = (PetscScalar*)malloc(sizeof(PetscScalar)*na*nb);

  PetscInt *idxs;
  idxs = (PetscInt*)malloc(sizeof(PetscInt)*na*nb);

  PetscInt idx = 0;
  for(int j = 0; j < nb; j++) 
    for(int i = 0; i < na; i++) {
      cs[i + na*j] = as[i] * bs[j] * c;
      idxs[i + na*j] = idx;
      idx++;
    }

  VecSetValues(*C, na*nb, idxs, cs, mode);

  VecRestoreArray(A, &as); VecRestoreArray(B, &bs); 
  free(cs); free(idxs);

  return 0;
}

PetscErrorCode VecSetSynthesize(Vec A, Vec B, PetscScalar c, 
				MPI_Comm comm, Vec *C){
  PetscErrorCode ierr;
  ierr = VecInitSynthesize(A, B, comm, C);
  ierr = VecSynthesize(A, B, c, C, INSERT_VALUES);
  VecAssemblyBegin(*C); VecAssemblyEnd(*C);
  return 0;
}

PetscErrorCode MatInitSynthesize(Mat A, Mat B, MPI_Comm comm, Mat *C) {

  PetscInt na, nb, ma, mb;
  PetscErrorCode ierr;
  ierr = MatGetSize(A, &na, &ma); CHKERRQ(ierr);
  ierr = MatGetSize(B, &nb, &mb); CHKERRQ(ierr);

  ierr = MatCreate(comm, C); CHKERRQ(ierr);
  ierr = MatSetSizes(*C, PETSC_DECIDE, PETSC_DECIDE, na*nb, ma*mb); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*C); CHKERRQ(ierr);
  ierr = MatSetUp(*C); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode MatSynthesize(Mat A, Mat B, PetscScalar c,
			     Mat *C, InsertMode mode) {
  PetscErrorCode ierr;
  PetscInt na, nb, ma, mb;
  ierr = MatGetSize(A, &na, &ma); CHKERRQ(ierr);
  ierr = MatGetSize(B, &nb, &mb); CHKERRQ(ierr);

  const PetscScalar **row_a, **row_b;
  PetscInt *ncols_a, *ncols_b;
  const PetscInt **cols_a, **cols_b;
  
  row_a = (const PetscScalar**)malloc(sizeof(PetscScalar*)*na);
  ncols_a = (PetscInt*)malloc(sizeof(PetscInt)*na);
  cols_a = (const PetscInt**)malloc(sizeof(PetscInt*)*na);
  row_b = (const PetscScalar**)malloc(sizeof(PetscScalar*)*nb);
  ncols_b = (PetscInt*)malloc(sizeof(PetscInt)*nb);
  cols_b = (const PetscInt**)malloc(sizeof(PetscInt*)*nb);

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
	  PetscScalar v = row_a[i_a][idx_a] * row_b[i_b][idx_b] * c;
	  ierr = MatSetValue(*C, i, j, v, mode);  CHKERRQ(ierr);
	}
      }
    }
  }

  for(int i = 0; i < na; i++)
    ierr = MatRestoreRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
  for(int i = 0; i < nb; i++)
    ierr = MatRestoreRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode MatSetSynthesize(Mat A, Mat B, PetscScalar c, 
				MPI_Comm comm, Mat *C) {
  PetscErrorCode ierr;
  ierr = MatInitSynthesize(A, B, comm, C); CHKERRQ(ierr);
  ierr = MatSynthesize(A, B, c, C, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*C, MAT_FINAL_ASSEMBLY);
  return 0;
}

PetscErrorCode MatInitSynthesize3(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D) {
  PetscInt na, nb, nc, ma, mb, mc;
  PetscErrorCode ierr;
  ierr = MatGetSize(A, &na, &ma); CHKERRQ(ierr);
  ierr = MatGetSize(B, &nb, &mb); CHKERRQ(ierr);
  ierr = MatGetSize(C, &nc, &mc); CHKERRQ(ierr);

  ierr = MatCreate(comm, D); CHKERRQ(ierr);
  ierr = MatSetSizes(*D, PETSC_DECIDE, PETSC_DECIDE, na*nb*nc, ma*mb*mc);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(*D); CHKERRQ(ierr);
  ierr = MatSetUp(*D); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode MatSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, 
			      Mat *D, InsertMode mode) {
  PetscErrorCode ierr;
  Mat BC;
  ierr = MatSetSynthesize(B, C, d, PETSC_COMM_SELF, &BC); CHKERRQ(ierr);
  
  ierr = MatSynthesize(A, BC, 1.0, D, mode); CHKERRQ(ierr);

  ierr = MatDestroy(&BC); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MatSetSynthesize3(Mat A, Mat B, Mat C, PetscScalar d, MPI_Comm comm, Mat *D) {
  PetscErrorCode ierr;
  ierr = MatInitSynthesize3(A, B, C, comm, D); CHKERRQ(ierr);
  ierr = MatSynthesize3(A, B, C, d, D, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*D, MAT_FINAL_ASSEMBLY);
  return 0;
}
