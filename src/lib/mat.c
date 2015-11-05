#include <rescol/mat.h>

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

  #if defined(PETSC_USE_COMPLEX)
  PetscReal a, b;
  while(fscanf(fp, "%d %d %lf %lf", &row, &col, &a, &b) != EOF) {
    dat = a + b * PETSC_i;
    ierr = MatSetValue(*mat, row, col, dat, INSERT_VALUES); CHKERRQ(ierr);
  }
  #else
  while(fscanf(fp, "%d %d %lf", &row, &col, &dat) != EOF) {
    ierr = MatSetValue(*mat, row, col, dat, INSERT_VALUES); CHKERRQ(ierr);
  }
  #endif

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
  PetscScalar eig, im_eig;
  PetscReal tol, error;
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
PetscErrorCode EPSCreateForBoundState(EPS *eps, MPI_Comm comm, Mat H, Mat S, PetscScalar target, EPSProblemType type) {

  EPSCreate(comm, eps);
  EPSSetOperators(*eps, H, S);
  EPSSetProblemType(*eps, type);
  EPSSetWhichEigenpairs(*eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(*eps, target);
  EPSSetFromOptions(*eps);

  return 0;
}


PetscErrorCode VecInitSynthesize(Vec A, Vec B, MPI_Comm comm, Vec *C) {

  PetscErrorCode ierr;
  PetscInt na, nb;
  ierr = VecGetSize(A, &na);  CHKERRQ(ierr);
  ierr = VecGetSize(B, &nb);  CHKERRQ(ierr);

  ierr =VecCreate(comm, C);  CHKERRQ(ierr);
  ierr =VecSetSizes(*C, PETSC_DECIDE, na*nb);  CHKERRQ(ierr);
  ierr =VecSetFromOptions(*C);  CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode VecSynthesize(Vec A, Vec B, PetscScalar c, 
			     Vec *C, InsertMode mode) {
  
  PetscInt na, nb;
  VecGetSize(A, &na); VecGetSize(B, &nb);

  PetscScalar *as, *bs, *cs;
  VecGetArray(A, &as); VecGetArray(B, &bs);
  PetscMalloc1(na*nb, &cs);

  PetscInt *idxs;
  PetscMalloc1(na*nb, &idxs);

  PetscInt idx = 0;
  for(int j = 0; j < nb; j++) 
    for(int i = 0; i < na; i++) {
      cs[i + na*j] = as[i] * bs[j] * c;
      idxs[i + na*j] = idx;
      idx++;
    }

  VecSetValues(*C, na*nb, idxs, cs, mode);

  VecRestoreArray(A, &as); VecRestoreArray(B, &bs); 
  PetscFree(idxs);
  PetscFree(cs);

  return 0;
}
PetscErrorCode VecSetSynthesize(Vec A, Vec B, PetscScalar c, 
				MPI_Comm comm, Vec *C){
  PetscErrorCode ierr;
  ierr = VecInitSynthesize(A, B, comm, C); CHKERRQ(ierr);
  ierr = VecSynthesize(A, B, c, C, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(*C);  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*C); CHKERRQ(ierr);
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
PetscErrorCode MatSetSynthesizeSlow(Mat A, Mat B, PetscScalar c, 
				MPI_Comm comm, Mat *C) {
  PetscErrorCode ierr;
  ierr = MatInitSynthesize(A, B, comm, C); CHKERRQ(ierr);
  ierr = MatSynthesize(A, B, c, C, INSERT_VALUES); CHKERRQ(ierr);
  MatAssemblyBegin(*C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*C, MAT_FINAL_ASSEMBLY);
  return 0;
}
PetscErrorCode MatSetSynthesizeFast(Mat A, Mat B, MPI_Comm comm, Mat *C) {
 
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

  PetscInt num_a = 0;
  PetscInt num_b = 0;
  for(int i = 0; i < na; i++) {
    ierr = MatGetRow(A, i, &ncols_a[i], &cols_a[i], &row_a[i]); CHKERRQ(ierr);
    num_a += ncols_a[i];
  }
  for(int i = 0; i < nb; i++) {
    ierr = MatGetRow(B, i, &ncols_b[i], &cols_b[i], &row_b[i]); CHKERRQ(ierr);
    num_b += ncols_b[i];
  }
  
  PetscInt num_c = num_a*num_b;
  PetscScalar *val;
  PetscInt *row, *col;
  val = (PetscScalar*)malloc(sizeof(PetscScalar)*num_c);
  row = (PetscInt*)malloc(sizeof(PetscInt)*num_c);
  col = (PetscInt*)malloc(sizeof(PetscInt)*num_c);

  PetscInt idx = 0;
  for(int i_a = 0; i_a < na; i_a++) {
    for(int idx_a = 0; idx_a < ncols_a[i_a]; idx_a++) {
      int j_a = cols_a[i_a][idx_a];
      for(int i_b = 0; i_b < nb; i_b++) { 
      	for(int idx_b = 0; idx_b < ncols_b[i_b]; idx_b++) {
	  int j_b = cols_b[i_b][idx_b];
	  row[idx] = i_a + i_b * na;
	  col[idx] = j_a + j_b * ma;
	  val[idx] = row_a[i_a][idx_a] * row_b[i_b][idx_b];
	  idx++;
	}
      }
    }
  }
  ierr = MatCreateSeqAIJFromTriple(comm, na*nb, ma*mb, row, col, val, C, num_c, 0);
  CHKERRQ(ierr);
  return 0;
}
PetscErrorCode MatSetSynthesize(Mat A, Mat B, PetscScalar c, 
				MPI_Comm comm, Mat *C) {
  return MatSetSynthesizeSlow(A, B, c, comm, C);
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
PetscErrorCode MatSetSynthesize3Fast(Mat A, Mat B, Mat C, MPI_Comm comm, Mat *D) {
  Mat BC;
  PetscErrorCode ierr;
  ierr = MatSetSynthesizeFast(B, C, comm, &BC); CHKERRQ(ierr);
  ierr = MatSetSynthesizeFast(A, BC, comm, D); CHKERRQ(ierr);
  MatDestroy(&BC);
  return 0;
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

PetscErrorCode LegGauss(int n, int i, PetscScalar* x, PetscScalar* w) {
  PetscScalar xs[45] = {
    0.0, -0.5773502691896257, 0.5773502691896257, -0.7745966692414834, 0.0, 0.7745966692414834, -0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526, -0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664, -0.932469514203152, -0.6612093864662645, -0.23861918608319693, 0.23861918608319693, 0.6612093864662645, 0.932469514203152, -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585, -0.9602898564975362, -0.7966664774136267, -0.525532409916329, -0.18343464249564984, 0.18343464249564984, 0.525532409916329, 0.7966664774136267, 0.9602898564975362, -0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
  PetscScalar ws[45] = {
    2.0, 1.0, 1.0, 0.5555555555555555, 0.888888888888889, 0.5555555555555555, 0.34785484513745396, 0.6521451548625462, 0.6521451548625462, 0.34785484513745396, 0.236926885056189, 0.4786286704993665, 0.5688888888888888, 0.4786286704993665, 0.236926885056189, 0.17132449237916944, 0.36076157304813883, 0.4679139345726918, 0.4679139345726918, 0.36076157304813883, 0.17132449237916944, 0.12948496616886826, 0.2797053914892774, 0.38183005050511937, 0.41795918367347007, 0.38183005050511937, 0.2797053914892774, 0.12948496616886826, 0.10122853629037527, 0.22238103445337512, 0.3137066458778874, 0.362683783378362, 0.362683783378362, 0.3137066458778874, 0.22238103445337512, 0.10122853629037527, 0.08127438836157427, 0.18064816069485748, 0.2606106964029357, 0.31234707704000275, 0.3302393550012596, 0.31234707704000275, 0.2606106964029357, 0.18064816069485748, 0.0812743883615742};

  *x = xs[n * (n-1)/2 + i];
  *w = ws[n * (n-1)/2 + i];
  return 0;
}
PetscErrorCode LobGauss(int n, int i, PetscScalar* x, PetscScalar* w) {
  PetscScalar xs[44] = {
    -1.0, 1.0, 
    -1.0, 0.0, 1.0, 
    -1.0, -0.4472135954999579, 0.4472135954999579, 1.0, 
    -1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0,
    -1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0, 
    -1.0, -0.830223896278567, -0.46884879347071423, 0.0, 0.46884879347071423,
 0.830223896278567, 1.0, 
    -1.0, -0.8717401485096066, -0.5917001814331423, -0.20929921790247888, 0.20929921790247888, 0.5917001814331423, 0.8717401485096066, 1.0,
    -1.0, -0.8997579954114602, -0.6771862795107377, -0.36311746382617816, 0.0,
 0.36311746382617816, 0.6771862795107377, 0.8997579954114602, 1.0};

  PetscScalar ws[44] = {
    1.0, 1.0, 
    0.3333333333333333, 1.3333333333333333, 0.3333333333333333, 
    0.16666666666666666, 0.8333333333333333, 0.8333333333333333, 0.16666666666666666, 0.1, 0.5444444444444444, 0.7111111111111111, 0.5444444444444444, 0.1, 0.06666666666666667, 0.3784749562978472, 0.5548583770354865, 0.5548583770354865, 0.3784749562978472, 0.06666666666666667, 0.047619047619047616, 0.2768260473615648, 0.4317453812098626, 0.4876190476190476, 0.4317453812098626, 0.2768260473615648, 0.047619047619047616, 0.03571428571428571, 0.21070422714350345, 0.3411226924835027, 0.4124587946587037, 0.4124587946587037, 0.3411226924835027, 0.21070422714350345, 0.03571428571428571, 0.027777777777777776, 0.1654953615608063, 0.2745387125001616, 0.34642851097304617, 0.37151927437641724, 0.34642851097304617, 0.2745387125001616, 0.1654953615608063, 0.027777777777777776};
  
  int idx = n*(n-1)/2 + i - 1;
  if(idx >= 44)
    SETERRQ(PETSC_COMM_SELF, 1, "index over flow");

  *x = xs[idx];
  *w = ws[idx];
  return 0;
}
