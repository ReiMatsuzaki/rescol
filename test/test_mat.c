#include <slepceps.h>
#include "unittest.h"
#include "../include/mat.h"
#include "../include/synthesize.h"

static char help[] = "Unit test for mat.c and synthesize.c \n\n";

PetscErrorCode testVecSplit() {
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;

  PetscScalar *x0s; PetscMalloc1(12, &x0s);
  for(int i = 0 ;i < 12; i++)
    x0s[i] = 1.1*i;
  Vec x; 
  VecCreateSeqWithArray(comm, 12, 12, x0s, &x);

  Vec *xs;
  ierr = VecGetSplit(x, 3, &xs); CHKERRQ(ierr);

  PetscScalar vs[4];
  PetscInt    idx[4] = {0, 1, 2, 3};
  ierr = VecGetValues(xs[0], 4, idx, vs); CHKERRQ(ierr);
  //  VecView(x, PETSC_VIEWER_STDOUT_SELF);
  //  VecView(xs[0], PETSC_VIEWER_STDOUT_SELF);
  ASSERT_DOUBLE_EQ(0.0, vs[0]);
  ASSERT_DOUBLE_EQ(1.1, vs[1]);
  ASSERT_DOUBLE_EQ(2.2, vs[2]);
  ASSERT_DOUBLE_EQ(3.3, vs[3]);

  ierr = VecGetValues(xs[1], 4, idx, vs); CHKERRQ(ierr);
  ASSERT_DOUBLE_EQ(4.4, vs[0]);
  ASSERT_DOUBLE_EQ(5.5, vs[1]);
  ASSERT_DOUBLE_EQ(6.6, vs[2]);
  ASSERT_DOUBLE_EQ(7.7, vs[3]);

  ierr = VecGetValues(xs[2], 4, idx, vs); CHKERRQ(ierr);
  ASSERT_DOUBLE_EQ(8.8, vs[0]);
  ASSERT_DOUBLE_EQ(9.9, vs[1]);
  ASSERT_DOUBLE_EQ(11.0, vs[2]);
  ASSERT_DOUBLE_EQ(12.1, vs[3]);
  
  VecRestoreSplit(x, 3, &xs);
  PetscFree(x0s);
  VecDestroy(&x);
  

  return 0;
}
int testMatMatDecomposedMult() {
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  Mat A, B;
  Vec x, y0, y1;
  MatCreate(comm, &A); MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2, 3); MatSetUp(A);
  MatCreate(comm, &B); MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, 4, 5); MatSetUp(B);
  VecCreate(comm, &x); VecSetSizes(x, PETSC_DECIDE, 15); VecSetUp(x);
  VecCreate(comm, &y0); VecSetSizes(y0, PETSC_DECIDE, 8); VecSetUp(y0);
  VecCreate(comm, &y1); VecSetSizes(y1, PETSC_DECIDE, 8); VecSetUp(y1);

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

  for(int i = 0; i < 15; i++)
    VecSetValue(x, i, 0.1+1.0*i, INSERT_VALUES);
  VecAssemblyBegin(x); VecAssemblyEnd(x);
  
  Mat C; 
  ierr = MatMatSynthesize(A, B, 1.2, MAT_INITIAL_MATRIX, &C); CHKERRQ(ierr);

  ierr = MatMult(C, x, y0); CHKERRQ(ierr);
  ierr = MatMatDecomposedMult(A, B, 1.2, x, y1); CHKERRQ(ierr);

  VecAXPY(y0, -1.0, y1);
  PetscReal norm;
  VecNorm(y0, NORM_1, &norm);
  ASSERT_DOUBLE_EQ(0.0, PetscRealPart(norm));
  
  MatDestroy(&A); MatDestroy(&B); MatDestroy(&C);
  VecDestroy(&x); VecDestroy(&y0); VecDestroy(&y1);

  return 0;
}
int testMatMatMatDecomposedMult() {
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscErrorCode ierr;
  Mat A; MatCreate(comm, &A); MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2, 3);
  Mat B; MatCreate(comm, &B); MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, 3, 4);
  Mat C; MatCreate(comm, &C); MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, 2, 5);
  Vec x; VecCreate(comm, &x); VecSetSizes(x, PETSC_DECIDE, 60);
  Vec y0;VecCreate(comm, &y0);VecSetSizes(y0,PETSC_DECIDE, 12);
  Vec y1;VecCreate(comm, &y1);VecSetSizes(y1,PETSC_DECIDE, 12);

  MatSetUp(A); MatSetUp(B); MatSetUp(C); VecSetUp(x); VecSetUp(y0); VecSetUp(y1);
  ierr = MatSetRandom(A, NULL); CHKERRQ(ierr);
  ierr = MatSetRandom(B, NULL); CHKERRQ(ierr);
  ierr = MatSetRandom(C, NULL); CHKERRQ(ierr);
  ierr = VecSetRandom(x, NULL); CHKERRQ(ierr);

  ierr = MatMatMatDecomposedMult(A, B, C, 1.2, x, y0);

  Mat D; MatMatMatSynthesize(A, B, C, 1.2, MAT_INITIAL_MATRIX, &D);
  ierr = MatMult(D, x, y1); CHKERRQ(ierr);

  VecAXPY(y0, -1.0, y1);
  PetscReal norm;
  VecNorm(y0, NORM_1, &norm);
  ASSERT_DOUBLE_NEAR(0.0, PetscRealPart(norm), pow(10.0, -10.0));

  MatDestroy(&A); MatDestroy(&B); MatDestroy(&C); MatDestroy(&D);
  VecDestroy(&x); VecDestroy(&y0); VecDestroy(&y1);
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
int testVecArrayLoad() {

/*
  MPI_Comm comm = PETSC_COMM_SELF;
  char path[256] = "tmpvec.dat";
  PetscErrorCode ierr;
  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(comm, path, FILE_MODE_READ, &viewer); CHKERRQ(ierr);  
  
  Vec *xs; 
  int n;
  VecArrayLoad(PETSC_VIEWER_STDOUT_SELF, &n, &xs);

*/  
  return 0;
}

int main(int argc, char **args) {
  
  SlepcInitialize(&argc, &args, (char*)0, help);
  MPI_Comm comm = PETSC_COMM_SELF;
  PrintTimeStamp(comm, "MatMatDecomposedMult", NULL);
  testMatMatDecomposedMult();
  PrintTimeStamp(comm, "MatMatMatDecomposedMult", NULL);
  testMatMatMatDecomposedMult();  
  PrintTimeStamp(comm, "VecSplit", NULL);
  testVecSplit();
  PrintTimeStamp(comm, "LegGauss", NULL);
  testLegGauss();
  PrintTimeStamp(comm, "LobGauss", NULL);
  testLobGauss();
  PrintTimeStamp(comm, "Coulomb", NULL);
  testPartialCoulomb();

  SlepcFinalize();
  return 0;
}
