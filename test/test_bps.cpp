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

