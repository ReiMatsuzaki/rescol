#include <slepceps.h>
#include "unittest.h"
#include <rescol/mat.h>


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
PetscErrorCode testMatSynthesize() {
  Mat A, B, C0, C1;
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

  MatSetSynthesizeFast(A, B, PETSC_COMM_WORLD, &C0);
  MatSetSynthesizeSlow(A, B, 1.0, PETSC_COMM_WORLD, &C1);

  const PetscScalar *row0, *row1;
  PetscInt ncols0, ncols1;
  const PetscInt *cols0, *cols1;
  MatGetRow(C0, 0, &ncols0, &cols0, &row0);
  MatGetRow(C1, 0, &ncols1, &cols1, &row1);
  ASSERT_EQ(ncols0, ncols1);
  ASSERT_EQ(cols0[0], cols1[0]);
  ASSERT_EQ(cols0[1], cols1[1]);

  return 0;
}
PetscErrorCode testMatSynthesize3() {
  Mat A, B, C0, C1;
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

  MatSetSynthesize3Fast(A, B, A, PETSC_COMM_WORLD, &C0);
  MatSetSynthesize3(A, B, A, 1.0, PETSC_COMM_WORLD, &C1);

  const PetscScalar *row0, *row1;
  PetscInt ncols0, ncols1;
  const PetscInt *cols0, *cols1;
  MatGetRow(C0, 0, &ncols0, &cols0, &row0);
  MatGetRow(C1, 0, &ncols1, &cols1, &row1);
  ASSERT_EQ(ncols0, ncols1);
  ASSERT_EQ(cols0[0], cols1[0]);
  ASSERT_EQ(cols0[1], cols1[1]);

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
int testLegGauss() {
  PetscScalar x, w;
  LegGauss(1, 0, &x, &w);
  ASSERT_DOUBLE_EQ(0.0, x);
  ASSERT_DOUBLE_EQ(2.0, w);

  LegGauss(2, 0, &x, &w);
  ASSERT_DOUBLE_EQ(-sqrt(1.0/3.0), x);
  ASSERT_DOUBLE_EQ(1.0, w);

  LegGauss(2, 1, &x, &w);
  ASSERT_DOUBLE_EQ(sqrt(1.0/3.0), x);
  ASSERT_DOUBLE_EQ(1.0, w);

  LegGauss(3, 0, &x, &w);
  ASSERT_DOUBLE_EQ(-sqrt(3.0/5.0), x);
  ASSERT_DOUBLE_EQ(5.0/9.0, w);

  LegGauss(3, 1, &x, &w);
  ASSERT_DOUBLE_EQ(0.0, x);
  ASSERT_DOUBLE_EQ(8.0/9.0, w);

  LegGauss(3, 2, &x, &w);
  ASSERT_DOUBLE_EQ(sqrt(3.0/5.0), x);
  ASSERT_DOUBLE_EQ(5.0/9.0, w);

  return 0;
}
int testLobGauss() {
  /* 
    (array([-1.        , -0.65465367,  0.        ,  0.65465367,  1.        ]),
    array([ 0.1       ,  0.54444444,  0.71111111,  0.54444444,  0.1       ]))
  */
  PetscScalar x, w; double e = pow(10.0, -8.0);
  LobGauss(2, 0, &x, &w);
  ASSERT_DOUBLE_EQ(-1.0, x); ASSERT_DOUBLE_EQ(1.0, w);
  LobGauss(2, 1, &x, &w);
  ASSERT_DOUBLE_EQ(1.0, x); ASSERT_DOUBLE_EQ(1.0, w);

  LobGauss(5, 0, &x, &w);
  ASSERT_DOUBLE_EQ(-1.0, x); ASSERT_DOUBLE_EQ(0.1, w);

  LobGauss(5, 1, &x, &w);
  ASSERT_DOUBLE_NEAR(-0.65465367, x, e); ASSERT_DOUBLE_NEAR(0.544444444444, w, e);

  LobGauss(5, 2, &x, &w);
  ASSERT_DOUBLE_EQ(0.0, x); ASSERT_DOUBLE_EQ(0.711111111111111, w);

  LobGauss(5, 3, &x, &w);
  ASSERT_DOUBLE_NEAR(0.65465367, x, e); ASSERT_DOUBLE_NEAR(0.544444444444, w, e);

  LobGauss(5, 4, &x, &w);
  ASSERT_DOUBLE_EQ(1.0, x); ASSERT_DOUBLE_EQ(0.1, w);
  
  return 0;
}
int testPartialCoulomb() {

  double v; 
  PartialCoulomb(0, 0.0, 1.1, &v);
  ASSERT_DOUBLE_EQ(1.0/1.1, v);

  PartialCoulomb(0, 1.1, 0.0, &v);
  ASSERT_DOUBLE_EQ(1.0/1.1, v);

  PartialCoulomb(1, 1.1, 0.0, &v);
  ASSERT_DOUBLE_EQ(0.0, v);
  PartialCoulomb(2, 1.1, 0.0, &v);
  ASSERT_DOUBLE_EQ(0.0, v);
  return 0;
}

int main(int argc, char **args) {
  
  SlepcInitialize(&argc, &args, (char*)0, help);

  PetscErrorCode ierr;
  ierr = testMat(); CHKERRQ(ierr);
  ierr = testMatSynthesize(); CHKERRQ(ierr);
  ierr = testMatSynthesize3(); CHKERRQ(ierr);  
  ierr = testVecSynthesize(); CHKERRQ(ierr);
  testLegGauss();
  testLobGauss();
  testPartialCoulomb();

  SlepcFinalize();
  return 0;
}
