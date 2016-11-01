#include <slepceps.h>
#include <gsl/gsl_sf_coupling.h>
#include <gtest/gtest.h>
#include "../include/angmoment.h"
#include "../include/y1s.h"
#include "../include/y2s.h"

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

  int jmax = 5;
  for(int j1 = 0; j1 < jmax; j1++)
    for(int j2 = 0; j2 < jmax; j2++) {
      int maxm = j1<j2 ? j1 : j2;
      for(int m = 0; m <= maxm; m++) {
	for(int q = 0; q <= jmax; q++) {
	  ASSERT_NEAR(Y1ElePq(j1, q, j2,
			      m,     m),
		      Y1ElePq(j2, q, j1,
			      m,     m),
		      0.000000001);
	}
      }
    }
}
TEST(Y1, Y1s_gerade) {
  
  Y1s y1s;
  Y1sCreate(PETSC_COMM_SELF, &y1s);
  Y1sSet(y1s, 0, GERADE, 5);
  int n; Y1sGetSize(y1s, &n);
  ASSERT_EQ(3, n);  
  if(getenv("SHOW_DEBUG"))
    Y1sView(y1s, PETSC_VIEWER_STDOUT_SELF);
  EXPECT_EQ(3, y1s->num);
  EXPECT_EQ(0, y1s->ls[0]);
  EXPECT_EQ(2, y1s->ls[1]);
  EXPECT_EQ(4, y1s->ls[2]);
  EXPECT_EQ(0, y1s->m);
  Y1sDestroy(&y1s);
  
}
TEST(Y1, Y1s_ungerade) {
  
  Y1s y1s;
  Y1sCreate(PETSC_COMM_SELF, &y1s);
  Y1sSet(y1s, 0, UNGERADE, 5);
  int n; Y1sGetSize(y1s, &n);
  ASSERT_EQ(3, n);
  if(getenv("SHOW_DEBUG"))
    Y1sView(y1s, PETSC_VIEWER_STDOUT_SELF);
  EXPECT_EQ(3, y1s->num);
  EXPECT_EQ(1, y1s->ls[0]);
  EXPECT_EQ(3, y1s->ls[1]);
  EXPECT_EQ(5, y1s->ls[2]);
  EXPECT_EQ(0, y1s->m);
  Y1sDestroy(&y1s);
  
}
TEST(Y1, Y1s_Y1s) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s a; Y1sCreate(comm, &a); Y1sSet(a, 0, GERADE, 6);
  Y1s b; Y1sCreate(comm, &b); Y1sSet(b, 0, UNGERADE, 5);

  PetscBool non0;
  Mat m;
  ierr = Y1sCreateY1MatOther(a, b, &m); ASSERT_EQ(0,ierr);
  ierr = Y1sYqkY1MatOther(a, b, 1, 0, m, &non0); ASSERT_EQ(0,ierr);

  ASSERT_TRUE(non0);

  int idx_i[1] = {0};
  int idx_j[1] = {1};
  PetscScalar xs[1];
  MatGetValues(m, 1, idx_i, 1, idx_j, xs);
  EXPECT_DOUBLE_EQ(0.0, PetscRealPart(xs[0]));

  Y1sDestroy(&a);
  Y1sDestroy(&b);
  MatDestroy(&m);
}
TEST(Y1, SetOne) {
  
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s y;
  ierr = Y1sCreate(comm, &y); ASSERT_EQ(0, ierr);
  ierr = Y1sSetOne(y, 3, 2); ASSERT_EQ(0, ierr);

  ASSERT_EQ(1, y->num);
  ASSERT_EQ(3, y->ls[0]);
  ASSERT_EQ(2, y->m);

}
/*
TEST(Y1, Pq_mat) {
  PetscErrorCode ierr;
  PetscBool non0;
  Y1s y1s;
  Mat pq_y; 
  ierr = Y1sCreateY1Mat(self->y1s, &pq_y); CHKERRQ(ierr);
  ierr = Y1sPqY1Mat(self->y1s, q, pq_y, &non0); CHKERRQ(ierr);
  MatAssemblyBegin(pq_y, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pq_y,   MAT_FINAL_ASSEMBLY);  
}
*/
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
  Y2sCreate(PETSC_COMM_SELF, &y2s);
  Y2sSet(y2s, SIGMA, GERADE, PLUS, 2);
  int n; Y2sGetSize(y2s, &n);
  ASSERT_EQ(5, n);
  Y2sDestroy(&y2s);
}

class TestY2s : public ::testing::Test {
public:
  Y2s y2s;
  
  virtual void SetUp() {
    Y2sCreate(PETSC_COMM_SELF, &y2s);
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

  if(getenv("SHOW_DEBUG"))
    Y2sView(y2s, PETSC_VIEWER_STDOUT_SELF);
}
TEST_F(TestY2s, S) {
  Mat M; Y2sCreateY2Mat(y2s, &M);
  Y2sSY2Mat(y2s, M);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(0, cols[0]);
  ASSERT_DOUBLE_EQ(1.0, PetscRealPart(row[0]));
  MatDestroy(&M);
}
TEST_F(TestY2s, Guess) {
  Vec v; Y2sCreateY2Vec(y2s, &v);
  Y2sGuessY2Vec(y2s, 0, 0, v);
  
  PetscScalar y[1];
  const PetscInt ix[1] = {0};
  VecGetValues(v, 1, ix, y);
  ASSERT_DOUBLE_EQ(1.0, PetscRealPart(y[0]));
  
  VecDestroy(&v);
}
TEST_F(TestY2s, Pq12) {
  Mat M; Y2sCreateY2Mat(y2s, &M);
  PetscBool find;
  Y2sPq12Y2Mat(y2s, 1, M, &find);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(1, cols[0]);
  ASSERT_DOUBLE_EQ(-0.57735026918962562, PetscRealPart(row[0]));
  MatDestroy(&M);
}
TEST_F(TestY2s, Pq1A) {
  Mat M; Y2sCreateY2Mat(y2s, &M);
  PetscBool find;
  Y2sPq1Y2Mat(y2s, 2, M, &find);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(4, cols[0]);
  ASSERT_DOUBLE_EQ(0.44721359549995793, PetscRealPart(row[0]));
  MatDestroy(&M);  

  Mat MM; Y2sCreateY2Mat(y2s, &MM);
  Y2sPq1Y2Mat(y2s, 1, MM, &find);
  ASSERT_FALSE(find);
}
TEST_F(TestY2s, Pq2A) {
  Mat M; Y2sCreateY2Mat(y2s, &M);
  PetscBool find;
  Y2sPq2Y2Mat(y2s, 2, M, &find);
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(M, 0, &ncols, &cols, &row);
  ASSERT_EQ(1, ncols);
  ASSERT_EQ(2, cols[0]);
  ASSERT_DOUBLE_EQ(0.44721359549995793, PetscRealPart(row[0]));
  MatDestroy(&M);  
}

int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
  
}
