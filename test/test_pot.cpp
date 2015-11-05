#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/pot.h>

static char help[] = "unit test for pot.c";

TEST(TestPOT, Harmonic) {
  POT harm; POTHarmCreate(&harm, 2.5);

#if defined(SHOW_DEBUG)
  POTView(harm);
#endif

  PetscScalar y;
  POTCalc(harm, 0.2, &y);
  ASSERT_DOUBLE_EQ(2.5*0.5*0.2*0.2, PetscRealPart(y));
		   
}

int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}
