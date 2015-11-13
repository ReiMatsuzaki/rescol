#include <slepceps.h>
#include "unittest.h"
#include <rescol/mat.h>


static char help[] = "Unit test for angmoment.c \n\n";

PetscErrorCode testVecSplit() {
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;

  PetscScalar *x0s; PetscMalloc1(12, &x0s);
  for(int i = 0 ;i < 12; i++)
    x0s[i] = 1.1*i;
  Vec x; 
  VecCreateSeqWithArray(comm, 12, 12, x0s, &x);

  Vec *xs; PetscMalloc1(3, &xs);
  ierr = VecSplit(x, 3, xs); CHKERRQ(ierr);

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
  
  for(int i = 0; i < 3; i++)
    VecDestroy(&xs[i]);
  PetscFree(xs);
  PetscFree(x0s);
  VecDestroy(&x);
  

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
  MPI_Comm comm = PETSC_COMM_SELF;
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
