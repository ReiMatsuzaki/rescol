#include <gtest/gtest.h>
#include <slepceps.h>
#include <rescol/oce1.h>
#include <rescol/eeps.h>

static char help[] = "Unit test for oce1.c";

TEST(OCE1Test, HAtom) {
  MPI_Comm comm = PETSC_COMM_SELF;
  Y1s y1s;    Y1sCreate(comm, &y1s); Y1sSetOne(y1s, 0, 1);
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
  int order = 5;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);
  OCE1 oce;   OCE1Create(comm, &oce); OCE1Set(oce, fem, y1s);
  Pot pot; PotCreate(comm, &pot); PotSetCoulombNE(pot, 0, 0.0, -1.0);
  
  Mat S; OCE1SMat(oce, MAT_INITIAL_MATRIX, &S, NULL);
  Mat H; OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  Mat V; OCE1PotMat(oce, ROT_SCALAR, pot, MAT_INITIAL_MATRIX, &V);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);
  
  EEPS eps; EEPSCreate(comm, &eps); 
  EEPSSetOperators(eps, H, S);
  EEPSSetTarget(eps, -0.2);
  EEPSSolve(eps);

  PetscScalar k;
  EPSGetEigenpair(eps->eps, 0, &k, 0, 0, 0);
  ASSERT_NEAR(-0.125, PetscRealPart(k), 0.00001);

  OCE1Destroy(&oce);
  PFDestroy(&pot);
  MatDestroy(&S); MatDestroy(&H); MatDestroy(&V);
  EEPSDestroy(&eps);
}
int _main(int argc, char **args) {
 ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  _main(argc, args);
  SlepcFinalize();
  return 0;
}
