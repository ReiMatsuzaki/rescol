#include <slepceps.h>
#include <gsl/gsl_sf_coupling.h>
#include <gtest/gtest.h>
#include <rescol/angmoment.h>

static char help[] = "Unit test for angmoment.c \n\n";

TEST(Utils, 6j) {
  ASSERT_DOUBLE_EQ(-1.0/3.0,
		   wigner6j(1, 1, 1, 
			    1, 1, 0));

}
TEST(Y1, Yqk) {
  
  // see ipython note named angmoment.ipynb
  ASSERT_NEAR(0.282094791774,
	      Y1EleYqk(+1, +0, +1,
		       +0, +0, +0),
	      0.000000001);
  ASSERT_NEAR(0.241795535806,
		     Y1EleYqk(+2, +2, +4,
			      +0, +0, +0),
		     0.000000001);
  ASSERT_NEAR(-0.0901118757864,
		     Y1EleYqk(+2, +2, +4,
			      +2, +1, +1),
		     0.000000001);
}
TEST(Y1, Pq) {
  // see ipython note named angmoment.ipynb
  ASSERT_DOUBLE_EQ(1.0,
		   Y1ElePq(+1, +0, +1,
			   +0,    +0));
  ASSERT_NEAR(0.383325939,
		     Y1ElePq(+2, +2, +4,
			     +0,     +0),
		     0.000000001);
  ASSERT_NEAR(0.349927106112,
		     Y1ElePq(+2, +2, +4,
			     +1,     +1),
		     0.000000001);
}
TEST(Y1, Y1s) {
  Y1s y1s;
  Y1sCreate(&y1s, PETSC_COMM_SELF);
  Y1sSet(y1s, SIGMA, GERADE, 5);
  int n; Y1sGetSize(y1s, &n);
  ASSERT_EQ(3, n);
  Y1sDestroy(&y1s);
}
TEST(Y2, Pq) {

  Y2 a = {1, 1, 1, 0};
  Y2 b = {1, 1, 1, 0};

  double eps = 0.000000001;
  ASSERT_NEAR(1.0, Y2ElePq12(a, 0, b), eps);
  ASSERT_NEAR(-0.2, Y2ElePq12(a, 2, b), eps);
  ASSERT_NEAR(1.0, Y2ElePq1A(a, 0, b), eps);
  ASSERT_NEAR(-0.2, Y2ElePq2A(a, 2, b), eps);

  Y2 c = {0, 0, 0, 0};
  Y2 d = {1, 1, 0, 0};
  ASSERT_NEAR(-0.57735026918962562, Y2ElePq12(c, 1, d), eps);
}
TEST(Y2, Set) {
  Y2s y2s;
  Y2sCreate(&y2s, PETSC_COMM_SELF);
  Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);
  int n; Y2sGetSize(y2s, &n);
  ASSERT_EQ(5, n);
  Y2sDestroy(&y2s);
}

class TestY2s : public ::testing::Test {
public:
  Y2s y2s;
  
  virtual void SetUp() {
    Y2sCreate(&y2s, PETSC_COMM_SELF);
    Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);
    /*
      (L1,L2,L, M)
      (0, 0, 0, 0)
      (1, 1, 0, 0)
      (0, 2, 2, 0)
      (1, 1, 2, 0)
      (2, 0, 2, 0)
     */
  }

  virtual void TearDown() {
    Y2sDestroy(&y2s);
  }

};
TEST_F(TestY2s, LMax) {
  int LMax;
  Y2sGetMaxL(y2s, &LMax);
  ASSERT_EQ(2, LMax);
}
TEST_F(TestY2s, S) {
  Mat M;
  Y2sSetSY2Mat(y2s, &M);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(0, cols[0]);
  ASSERT_DOUBLE_EQ(1.0, row[0]);
  MatDestroy(&M);
}
TEST_F(TestY2s, Guess) {
  Vec v;
  Y2sSetGuessY2Vec(y2s, 0, 0, &v);
  
  PetscScalar y[1];
  const PetscInt ix[1] = {0};
  VecGetValues(v, 1, ix, y);
  ASSERT_DOUBLE_EQ(1.0, y[0]);
  
  VecDestroy(&v);
}
TEST_F(TestY2s, Pq12) {
  Mat M;
  Y2sSetPq12Y2Mat(y2s, 1, &M);
  //  Y2sView(y2s);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(1, cols[0]);
  ASSERT_DOUBLE_EQ(-0.57735026918962562, row[0]);
  MatDestroy(&M);
}
TEST_F(TestY2s, Pq1A) {
  Mat M;
  Y2sSetPq1AY2Mat(y2s, 2, &M);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(4, cols[0]);
  ASSERT_DOUBLE_EQ(0.44721359549995793, row[0]);
  MatDestroy(&M);  

  Mat MM;
  Y2sSetPq1AY2Mat(y2s, 1, &MM);
  ASSERT_TRUE(MM==NULL);
}
TEST_F(TestY2s, Pq2A) {
  Mat M;
  Y2sSetPq2AY2Mat(y2s, 2, &M);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(2, cols[0]);
  ASSERT_DOUBLE_EQ(0.44721359549995793, row[0]);
  MatDestroy(&M);  
}

/*
PetscErrorCode testY2ElePq() {

  PetscErrorCode ierr;
  CoupledY a = {1, 1, 1, 0};
  CoupledY b = {1, 1, 1, 0};
  ASSERT_DOUBLE_EQ(1.0, Y2ElePqR12(a, 0, b));
  ASSERT_DOUBLE_EQ(-0.2, Y2ElePqR12(a, 2, b));
  ASSERT_DOUBLE_EQ(-0.2, Y2ElePqR1A(a, 0, b));
  ASSERT_DOUBLE_EQ(-0.2, Y2ElePqR2A(a, 2, b));
  return 0;

}


petscErrorCode testY2Set() {

  PetscErrorCode ierr;
  Y2 y2s;
  ierr = Y2SetCreate(&y2s, PETSC_COMM_SELF); CHKERRQ(ierr);
  ierr = Y2SetSet(y2s, 0, Y2GERADE, Y2PLUS, 2); CHKERRQ(ierr);
  int n;
  ierr = Y2SetGetSize(y2s, &n); CHKERRQ(ierr);
  ASSERT_EQ(5, n);
  ierr = Y2SetDestroy(&y2s); CHKERRQ(ierr);
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
*/
int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  
  /*
  testY1ElePq();
  testY1s();
  testY2ElePq();
  testY2s();
  */

  // ierr = testY2ElePq(); CHKERRQ(ierr);
  // ierr = testY1s_Lambda(); CHKERRQ(ierr);

  SlepcFinalize();
  return 0;
  
}
