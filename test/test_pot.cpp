#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/pot.h>

static char help[] = "unit test for pot.c";

TEST(TestPOT, Harmonic) {
  PF harm; PFCreate(MPI_COMM_SELF, 1, 1, &harm); PFSetHarm(harm, 2.5);

  if(getenv("SHOW_DEBUG"))
    PFView(harm, PETSC_VIEWER_STDOUT_SELF);

  PetscScalar y[1];
  PetscScalar x[1] = {0.2};
  PFApply(harm, 1, x, y);
  ASSERT_DOUBLE_EQ(2.5*0.5*0.2*0.2, PetscRealPart(y[0]));
		   
}

int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}
