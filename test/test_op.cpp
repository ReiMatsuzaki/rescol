#include <petscmat.h>
#include <gtest/gtest.h>
#include <rescol/bspline.h>
#include <rescol/pot.h>

static char help[] = "Unit test for r1op.c, ry1op.c \n\n";

TEST(first, fist) {

  ASSERT_EQ(2, 1+1);

}
TEST(R1Op, D2) {
  MPI_Comm comm = PETSC_COMM_SELF;
  R1Op d2; R1OpCreate(comm, &d2); 
  R1OpSet(d2, R1OpD2, NULL, NULL);
  
  BPS bps;BPSCreate(comm, &bps);BPSSetLine(bps, 10, 10.0);
  BSS bss;BSSCreate(comm, &bss);BSSSetKnots(bss, 4, bps);
  BSSSetUp(bss);

  Mat M1; BSSCreateMat(bss, 1, &M1);
  Mat M2; BSSCreateMat(bss, 1, &M2);

  BSSD2R1Mat(bss, M1);
  BSSOpMat(bss, d2, M2);

  MatAXPY(M1, -1.0, M2);
  PetscReal a;
  MatAbs(M1, a);
  ASSERT_DOUBLE_EQ(0.0, a);

}
TEST(R1Op, Vne) {
  MPI_Comm comm = PETSC_COMM_SELF;
  PF pf; PFCreate(comm, 1, 1, &pf); 
  PFSetCoulombNE(pf, 0, 0.0);
  R1Op vne; R1OpCreate(comm, &vne); 
  R1OpSet(vne, R1OpPF, vne, NULL);
  
  BPS bps;BPSCreate(comm, &bps);BPSSetLine(bps, 10, 10.0);
  BSS bss;BSSCreate(comm, &bss);BSSSetKnots(bss, 4, bps);
  BSSSetUp(bss);

  Mat M1; BSSCreateMat(bss, 1, &M1);
  Mat M2; BSSCreateMat(bss, 1, &M2);

  BSSD2R1Mat(bss, M1);
  BSSOpMat(bss, vne, M2);

  MatAXPY(M1, -1.0, M2);
  PetscReal a;
  MatAbs(M1, a);
  ASSERT_DOUBLE_EQ(0.0, a);
}

int main(int argc, char **args) {

  PetscInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  PetscFinalize();
  return 0;

}
