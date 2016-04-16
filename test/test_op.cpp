#include <petscmat.h>
#include <gtest/gtest.h>
#include <rescol/fem_inf.h>
#include <rescol/pot.h>
#include <rescol/op.h>

static char help[] = "Unit test for r1op.c, ry1op.c \n\n";

TEST(first, fist) {

  ASSERT_EQ(2, 1+1);

}

class TestOp : public ::testing::Test {
public:
  MPI_Comm comm;
  BSS bss;
  virtual void SetUp() {
    PetscErrorCode ierr;
    comm = PETSC_COMM_SELF;
    BPS bps; 
    ierr = BPSCreate(comm, &bps); ASSERT_EQ(0, ierr);
    ierr = BPSSetLine(bps, 10, 10.0); ASSERT_EQ(0, ierr);
    ierr = BSSCreate(comm, &bss);ASSERT_EQ(0, ierr);
    ierr = BSSSetKnots(bss, 4, bps);ASSERT_EQ(0, ierr);
    ierr = BSSSetUp(bss);ASSERT_EQ(0, ierr);
  }
  virtual void TearDown() {
    PetscErrorCode ierr;
    ierr =BSSDestroy(&bss);ASSERT_EQ(0, ierr);
  }
};
TEST_F(TestOp, D2) {
  PetscErrorCode ierr;

  Op d2;
  ierr = OpCreate(comm, &d2); ASSERT_EQ(0, ierr);
  ierr = OpSetD2(d2); ASSERT_EQ(0, ierr);

  if(getenv("SHOW_DEBUG")) {
    ierr = OpView(d2, PETSC_VIEWER_STDOUT_SELF); ASSERT_EQ(0, ierr);
  }

  Mat M1; 
  ierr = BSSCreateR1Mat(this->bss, &M1); ASSERT_EQ(0, ierr);
  ierr = BSSD2R1Mat(bss, M1); ASSERT_EQ(0, ierr);

  Mat M2; 
  ierr = BSSCreateR1Mat(bss, &M2); ASSERT_EQ(0, ierr);
  ierr = BSSOpMat(bss, d2, M2); ASSERT_EQ(0, ierr);

  MatAXPY(M1, -1.0, M2, DIFFERENT_NONZERO_PATTERN);
  PetscReal a;
  MatNorm(M1, NORM_1, &a);
  ASSERT_DOUBLE_EQ(0.0, a);

  OpDestroy(&d2);
  MatDestroy(&M1);
  MatDestroy(&M2);
}
TEST_F(TestOp, Vne) {
  PetscErrorCode ierr;
  PF pf;  PotCreate(comm, &pf); PotSetCoulombNE(pf, 1, 0.0, 1.0);
  Op vne; OpCreate(comm, &vne); OpSetPF(vne, pf);
  
  Mat M1; BSSCreateR1Mat(bss, &M1); BSSENR1Mat(bss, 1, 0.0, M1);
  Mat M2; BSSCreateR1Mat(bss, &M2); BSSOpMat(bss, vne, M2);

  if(getenv("SHOW_DEBUG")) {
    ierr = OpView(vne, PETSC_VIEWER_STDOUT_SELF); ASSERT_EQ(0, ierr);
  }

  ierr = MatAXPY(M1, -1.0, M2, DIFFERENT_NONZERO_PATTERN); ASSERT_EQ(0, ierr);
  PetscReal a;
  ierr = MatNorm(M1, NORM_1, &a);ASSERT_EQ(0, ierr);
  ASSERT_DOUBLE_EQ(0.0, a);

  MatDestroy(&M1); MatDestroy(&M2);
  OpDestroy(&vne);
}

int _main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
int main(int argc, char **args) {

  PetscInitialize(&argc, &args, (char*)0, help);
  _main(argc, args);
  PetscFinalize();
  return 0;

}
