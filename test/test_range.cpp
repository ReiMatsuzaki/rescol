#include <iostream>
#include <slepceps.h>
#include <gtest/gtest.h>
#include "../include/range.h"

using namespace std;

static char help[] = "Test for range.c\n\n";

TEST(TestRange, first) {

  PetscErrorCode ierr;
  Range range;
  ierr = RangeCreate(PETSC_COMM_SELF, &range);
  ASSERT_EQ(0, ierr);
  //  ierr = RangeSet(range, 1.0, 5.0, 5);
  ierr = RangeSetFromStr(range, "1.0:5.0:5");
  ASSERT_EQ(0, ierr);
  ierr = RangeSetName(range, "xs");
  ASSERT_EQ(0, ierr);

  PetscReal x_ref(1.0);
  PetscReal x = RangeInit(range);
  while(RangeNext(range, &x)) {
    ASSERT_DOUBLE_EQ(x_ref, x);
    x_ref += 1.0;
  }
  ierr = RangeView(range, PETSC_VIEWER_STDOUT_SELF);
  ASSERT_EQ(0, ierr);
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
