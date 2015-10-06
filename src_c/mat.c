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

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat) {
  FILE* fp;
  PetscInt i, col, row;
  PetscErrorCode ierr;
  int num_data, num_row, num_col;
  PetscScalar dat;


  if((fp = fopen(path, "r")) == NULL) {
    const char* msg = "Failed to open file.\0";
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, msg);
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
  MatCreate(PETSC_COMM_WORLD, mat);
  MatSetSizes(*mat, PETSC_DECIDE, PETSC_DECIDE, num_row, num_col);
  MatSetFromOptions(*mat);
  MatSetUp(*mat);

  i = 0;
  while(fscanf(fp, "%d %d %lf", &row, &col, &dat) != EOF) {
    i++;
    ierr = MatSetValue(*mat, row, col, dat, INSERT_VALUES); CHKERRQ(ierr);
  }
  fclose(fp);
  MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);
  return ierr;
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
