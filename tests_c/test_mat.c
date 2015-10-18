#include <slepceps.h>
#include "unittest.h"
#include "../src_c/mat.h"


static char help[] = "Unit test for angmoment.c \n\n";

PetscErrorCode testMat() {

  Mat A, B, C;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatCreate(PETSC_COMM_WORLD, &B);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2, 3);
  MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, 4, 5);
  MatSetFromOptions(A);
  MatSetFromOptions(B);
  MatSetUp(A);
  MatSetUp(B);

  MatSetValue(A, 0, 0, 1.0, INSERT_VALUES);
  MatSetValue(A, 1, 0, 3.0, INSERT_VALUES);
  MatSetValue(A, 1, 1, 1.0, INSERT_VALUES);
  MatSetValue(A, 1, 2, 2.0, INSERT_VALUES);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatSetValue(B, 0, 0, 1.0, INSERT_VALUES);
  MatSetValue(B, 1, 0, 3.0, INSERT_VALUES);
  MatSetValue(B, 1, 1, 1.0, INSERT_VALUES);
  MatSetValue(B, 1, 2, 2.0, INSERT_VALUES);
  MatSetValue(B, 3, 4, 1.0, INSERT_VALUES);
  MatSetValue(B, 0, 3, 2.0, INSERT_VALUES);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

  MatInitSynthesize(A, B, PETSC_COMM_WORLD, &C);
  MatSynthesize(A, B, 0.1, &C, INSERT_VALUES);
  MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(C, 0, &ncols, &cols, &row);
  ASSERT_EQ(ncols, 2);
  ASSERT_DOUBLE_EQ(0.1, row[0]);   ASSERT_DOUBLE_EQ(0.2, row[1]);
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(9, cols[1]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  MatGetRow(C, 1, &ncols, &cols, &row);
  ASSERT_EQ(6, ncols);
  ASSERT_DOUBLE_EQ(0.3, row[0]);   ASSERT_DOUBLE_EQ(0.1, row[1]);
  ASSERT_DOUBLE_EQ(0.2, row[2]);   ASSERT_DOUBLE_EQ(0.6, row[3]);
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(1, cols[1]);
  ASSERT_EQ(2, cols[2]);   ASSERT_EQ(9, cols[3]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  MatGetRow(C, 3, &ncols, &cols, &row);
  ASSERT_EQ(ncols, 9);
  ASSERT_DOUBLE_EQ(0.9, row[0]);   ASSERT_DOUBLE_EQ(0.3, row[1]);
  ASSERT_DOUBLE_EQ(0.6, row[2]);   ASSERT_DOUBLE_EQ(0.3, row[3]);
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(1, cols[1]);
  ASSERT_EQ(2, cols[2]);   ASSERT_EQ(3, cols[3]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  /*
    A = 1 0 0
        3 1 2
    B = 1 0 0 2 0 
        3 1 2 0 0
	0 0 0 0 0
	0 0 0 0 1 
    C = 100 000 000 200 000
        312 000 000 624 000

	300 100 200 000 000
	936 312 624 000 000
	....
   */

  MatDestroy(&A);
  MatDestroy(&B);
  
  return 0;
}
PetscErrorCode testVecSynthesize() {

  Vec A, B;
  MPI_Comm comm = PETSC_COMM_WORLD;
  VecCreate(comm, &A); VecCreate(comm, &B);
  VecSetSizes(A, PETSC_DECIDE, 3); VecSetSizes(B, PETSC_DECIDE, 2);
  VecSetFromOptions(A); VecSetFromOptions(B);

  VecSetValue(A, 0, 1.0, INSERT_VALUES);
  VecSetValue(A, 1, 2.0, INSERT_VALUES);
  VecSetValue(A, 2, 3.0, INSERT_VALUES);

  VecSetValue(B, 0, 0.3, INSERT_VALUES);
  VecSetValue(B, 1, 2.2, INSERT_VALUES);

  VecAssemblyBegin(A); VecAssemblyEnd(A);
  VecAssemblyBegin(B); VecAssemblyEnd(B);

  Vec C;
  VecSetSynthesize(A, B, 0.2, comm, &C);

  PetscScalar *cs;
  PetscInt n;
  VecGetSize(C, &n);
  VecGetArray(C, &cs);  
  
  ASSERT_EQ(6, n);
  ASSERT_DOUBLE_EQ(1.0*0.3*0.2, cs[0]);
  ASSERT_DOUBLE_EQ(2.0*0.3*0.2, cs[1]);
  ASSERT_DOUBLE_EQ(3.0*0.3*0.2, cs[2]);
  ASSERT_DOUBLE_EQ(1.0*2.2*0.2, cs[3]);
  ASSERT_DOUBLE_EQ(2.0*2.2*0.2, cs[4]);
  ASSERT_DOUBLE_EQ(3.0*2.2*0.2, cs[5]);

  VecRestoreArray(C, &cs);
  VecDestroy(&A); VecDestroy(&B);

  return 0;
}

int main(int argc, char **args) {
  
  SlepcInitialize(&argc, &args, (char*)0, help);

  PetscErrorCode ierr;
  ierr = testMat(); CHKERRQ(ierr);
  ierr = testVecSynthesize(); CHKERRQ(ierr);

  SlepcFinalize();
  return 0;
}
