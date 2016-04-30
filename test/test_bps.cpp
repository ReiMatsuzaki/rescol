#include <slepceps.h>
#include <stdlib.h>
#include <rescol/bps.h>
#include <gtest/gtest.h>

static char help[] = "Unit test for dvr.c \n\n";

TEST(TestBPS, Line) {

  PetscErrorCode ierr;
  BPS bps;
  ierr = BPSCreate(PETSC_COMM_SELF, &bps); ASSERT_EQ(0, ierr);
  ierr = BPSSetLine(bps, 5.0, 6); ASSERT_EQ(0, ierr);

  double *zs;
  int num;
  ierr = BPSGetZs(bps, &zs, &num);

  ASSERT_EQ(6, num);
  
  for(int i = 0; i < num; i++)
    ASSERT_DOUBLE_EQ(i*1.0, zs[i]);

  PetscFree(zs);
  BPSDestroy(&bps);
}
TEST(TestBPS, Exp) {

  PetscErrorCode ierr;
  BPS bps;
  ierr = BPSCreate(PETSC_COMM_SELF, &bps); ASSERT_EQ(0, ierr);
  ierr = BPSSetExp(bps, 3.3, 6, 5.0); ASSERT_EQ(0, ierr);

  double *zs;
  int num;
  ierr = BPSGetZs(bps, &zs, &num);

  ASSERT_EQ(6, num);

  if(getenv("SHOW_DEBUG"))
    BPSView(bps, PETSC_VIEWER_STDOUT_SELF);

  PetscFree(zs);
  BPSDestroy(&bps);
}
TEST(TestBPS, Getter) {

  PetscErrorCode ierr;
  BPS bps;
  ierr = BPSCreate(PETSC_COMM_SELF, &bps);  ASSERT_EQ(0, ierr);
  ierr = BPSSetLine(bps, 5.0, 6);  ASSERT_EQ(0, ierr);

  PetscReal z0, z1;
  ierr = BPSGetEdge(bps, 1, &z0, &z1); ASSERT_EQ(0, ierr);
  EXPECT_DOUBLE_EQ(1.0, z0);
  EXPECT_DOUBLE_EQ(2.0, z1);

  PetscBool is_in_q;
  BPSInElementQ(bps, 0, -0.1, &is_in_q); EXPECT_FALSE(is_in_q);
  BPSInElementQ(bps, 0, 0.0, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 0, 0.5, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 0, 1.0, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 0, 1.1, &is_in_q); EXPECT_FALSE(is_in_q);

  BPSInElementQ(bps, 2, 1.9, &is_in_q); EXPECT_FALSE(is_in_q);
  BPSInElementQ(bps, 2, 2.0, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 2, 2.4, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 2, 3.0, &is_in_q); EXPECT_TRUE(is_in_q);
  BPSInElementQ(bps, 2, 3.1, &is_in_q); EXPECT_FALSE(is_in_q);

  BPSDestroy(&bps);

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

