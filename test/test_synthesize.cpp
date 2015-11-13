#include <gtest/gtest.h>
#include <petscmat.h>
#include <rescol/synthesize.h>

static char help[] = "Unit test for synthesize.c";

TEST(First, first) {

  ASSERT_EQ(2, 1+1);

}
TEST(TestVecSynthesize, TwoVec) {
  MPI_Comm comm = PETSC_COMM_SELF;
  Vec a; VecCreate(comm, &a); VecSetSizes(a, PETSC_DECIDE, 3); VecSetUp(a);
  Vec b; VecCreate(comm, &b); VecSetSizes(b, PETSC_DECIDE, 2); VecSetUp(b);
  for(int i = 0; i < 3; i++)
    VecSetValue(a, i, (i+1)*1.0, INSERT_VALUES);
  for(int i = 0; i < 2; i++)
    VecSetValue(b, i, (i+1)*2.0, INSERT_VALUES);

  Vec c; VecVecSynthesize(a, b, 1.1, MAT_INITIAL_MATRIX, &c);
  PetscScalar *vs;
  VecGetArray(c, &vs);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[0]), 2.2);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[1]), 4.4);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[2]), 6.6);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[3]), 4.4);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[4]), 8.8);
  ASSERT_DOUBLE_EQ(PetscRealPart(vs[5]), 13.2);
  VecRestoreArray(c, &vs);
  VecDestroy(&a);VecDestroy(&b);VecDestroy(&c);

}
TEST(TestMatSynthesize, TwoMat) {
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

  //  MatInitSynthesize(A, B, PETSC_COMM_WORLD, &C);
  MatMatSynthesize(A, B, 0.1, MAT_INITIAL_MATRIX, &C);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(C, 0, &ncols, &cols, &row);
  ASSERT_EQ(ncols, 2);
  ASSERT_DOUBLE_EQ(0.1, PetscRealPart(row[0]));   
  ASSERT_DOUBLE_EQ(0.2, PetscRealPart(row[1]));
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(9, cols[1]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  MatGetRow(C, 1, &ncols, &cols, &row);
  ASSERT_EQ(6, ncols);
  ASSERT_DOUBLE_EQ(0.3, PetscRealPart(row[0])); 
  ASSERT_DOUBLE_EQ(0.1, PetscRealPart(row[1]));
  ASSERT_DOUBLE_EQ(0.2, PetscRealPart(row[2]));
  ASSERT_DOUBLE_EQ(0.6, PetscRealPart(row[3]));
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(1, cols[1]);
  ASSERT_EQ(2, cols[2]);   ASSERT_EQ(9, cols[3]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  MatGetRow(C, 3, &ncols, &cols, &row);
  ASSERT_EQ(ncols, 9);
  ASSERT_DOUBLE_EQ(0.9, PetscRealPart(row[0])); 
  ASSERT_DOUBLE_EQ(0.3, PetscRealPart(row[1]));
  ASSERT_DOUBLE_EQ(0.6, PetscRealPart(row[2]));
  ASSERT_DOUBLE_EQ(0.3, PetscRealPart(row[3]));
  ASSERT_EQ(0, cols[0]);   ASSERT_EQ(1, cols[1]);
  ASSERT_EQ(2, cols[2]);   ASSERT_EQ(3, cols[3]);
  MatRestoreRow(C, 0, &ncols, &cols, &row);

  MatDestroy(&A); MatDestroy(&B); MatDestroy(&C);   

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
}
TEST(TestMatSynthesize, ThreeMat) {
  Mat A, B, C0;
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

  MatMatMatSynthesize(A, B, A, 1.0, MAT_INITIAL_MATRIX, &C0);

  Mat AB; MatMatSynthesize(A, B, 1.0, MAT_INITIAL_MATRIX, &AB);
  Mat C1; MatMatSynthesize(AB,A, 1.0, MAT_INITIAL_MATRIX, &C1);

  const PetscScalar *row0, *row1;
  PetscInt ncols0, ncols1;
  const PetscInt *cols0, *cols1;
  MatGetRow(C0, 0, &ncols0, &cols0, &row0);
  MatGetRow(C1, 0, &ncols1, &cols1, &row1);
  ASSERT_EQ(ncols0, ncols1);
  ASSERT_EQ(cols0[0], cols1[0]);
  ASSERT_EQ(cols0[1], cols1[1]);

  MatDestroy(&A); MatDestroy(&B); MatDestroy(&C0); MatDestroy(&C1);
  MatDestroy(&AB);
}

int _main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
int main (int argc, char **args) {
  PetscInitialize(&argc, &args, (char*)0, help);
  _main(argc, args);
  PetscFinalize();
  return 0;
}

