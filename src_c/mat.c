#include "mat.h"

PetscErrorCode MatCreateFromCOOFormatFile(char* path, Mat* mat) {
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
