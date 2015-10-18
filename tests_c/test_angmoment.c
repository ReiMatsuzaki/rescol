#include <slepceps.h>
#include "unittest.h"
#include "../src_c/angmoment.h"
#include <Python/Python.h>

static char help[] = "Unit test for angmoment.c \n\n";

int testPyListToMat() {

  PyObject *pModule, *pFunc, *pArgs, *pList, *pValue;
  Py_Initialize();

  char path[2] = ".\0";
  PySys_SetPath(path);
  pModule = PyImport_ImportModule("for_test"); 
  
  ASSERT_NON_NULL(pModule);
  
  pFunc = PyObject_GetAttrString(pModule, "func1");
  pArgs = PyTuple_New(0);

  pList = PyObject_CallObject(pFunc, pArgs);
  Py_DECREF(pArgs); Py_DECREF(pFunc); Py_DECREF(pModule);

  ASSERT_EQ(3, (int)PyList_Size(pList));

  pValue = PyList_GetItem(pList, 0);
  ASSERT_DOUBLE_EQ(1.1, PyFloat_AsDouble(pValue));
  Py_DECREF(pValue);

  pValue = PyList_GetItem(pList, 1);
  ASSERT_DOUBLE_EQ(2.3, PyFloat_AsDouble(pValue));
  Py_DECREF(pValue);

  Py_DECREF(pList);
  return 0;
}
PetscErrorCode testY1Mat_Yqk() {
  
  ASSERT_DOUBLE_EQ(Y1Mat_Yqk(2, 1, 2,
			     -1, 0, 1),
		   PyGaunt(2, 1, 2,
			    -1, 0, 1));
  return 0;
}
PetscErrorCode testY1s_Pq() {
  Y1s ys;
  PetscErrorCode ierr;

  ierr = Y1sCreate(&ys, 0, 4, GERADE, 1); CHKERRQ(ierr);

  Mat A;
  Y1sCreateY1Mat(ys, PETSC_COMM_WORLD, &A);
  Y1sCalcPqY1Mat(ys, 1, A, INSERT_VALUES);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  ierr =   MatGetRow(A, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(1, cols[0]);
  ASSERT_DOUBLE_EQ(0.4886025119029199, row[0])
  ierr = MatRestoreRow(A, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr =   MatGetRow(A, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(2, cols[0]);
  ASSERT_DOUBLE_EQ(0.4886025119029199, row[0])
  ierr = MatRestoreRow(A, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  Mat B;
  Y1sCreateY1Mat(ys, PETSC_COMM_WORLD, &B);
  Y1sCalcPqY1Mat(ys, 1, B, INSERT_VALUES);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

  ierr =   MatGetRow(B, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(2, ncols);
  ASSERT_EQ(1, cols[0]); ASSERT_EQ(2, cols[1]);
  ASSERT_DOUBLE_EQ(-0.09011187578643429, row[0]);
  ASSERT_DOUBLE_EQ(-0.2207281154418226, row[1]);
  ierr = MatRestoreRow(B, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr =   MatGetRow(B, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(2, ncols);
  ASSERT_EQ(1, cols[0]); ASSERT_EQ(2, cols[1]);
  ASSERT_DOUBLE_EQ(-0.2207281154418226, row[0]);
  ASSERT_DOUBLE_EQ(-0.13926380803358027, row[1]);
  ierr = MatRestoreRow(B, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  //q == 1 =>
  //[[0.0, 0.0, 0.0], 
  // [0.0, 0.4886025119029199, 0.0], 
  // [0.0, 0.0, 0.4886025119029199]]

  // q == 2 =>
  // [[0.0, 0.0, 0.0], 
  // [0.0, -0.09011187578643429, -0.2207281154418226], 
  // [0.0, -0.2207281154418226, -0.13926380803358027]]  
  return 0;
}
PetscErrorCode testY1s_Lambda() {
  PetscErrorCode ierr;
  Y1s ys;
  // l0=0, L1=5, Gerade sym, M=0
  ierr = Y1sCreate(&ys, 0, 6, GERADE, 0); CHKERRQ(ierr);

  ASSERT_EQ(3, ys->num);
  ASSERT_EQ(0, ys->ls[0]);
  ASSERT_EQ(2, ys->ls[1]);
  ASSERT_EQ(4, ys->ls[2]);

  Mat M;
  Y1sCreateY1Mat(ys, PETSC_COMM_WORLD, &M);
  Y1sCalcLambdaY1Mat(ys, M, INSERT_VALUES);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

  int n, m;
  MatGetSize(M, &n, &m);
  ASSERT_EQ(3, n);
  ASSERT_EQ(3, m);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(0, ncols);
  MatRestoreRow(M, 0, &ncols, &cols, &row);

  MatGetRow(M, 1, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(1, cols[0]);
  ASSERT_DOUBLE_EQ(6.0, row[0]);
  MatRestoreRow(M, 0, &ncols, &cols, &row);
  
  MatDestroy(&M);
  Y1sDestroy(&ys);
  
  return 0;
}

int main(int argc, char **args) {
  
  SlepcInitialize(&argc, &args, (char*)0, help);

  printf("AABBCC\n");
  testPyListToMat();

  PetscErrorCode ierr;
  ierr = testY1Mat_Yqk(); CHKERRQ(ierr);
  ierr = testY1s_Pq(); CHKERRQ(ierr);
  //ierr = testY1s_Lambda(); CHKERRQ(ierr);
  SlepcFinalize();
  return 0;
}
