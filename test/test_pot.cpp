#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/pot.h>

static char help[] = "unit test for pot.c";

TEST(String, str_nele) {
  EXPECT_EQ(2, str_nele("aaa|bbbb|ccc", '|'));
}
void test_str2pot(const char str[], Pot p_ref, double x0, PetscErrorCode* p_ierr) {

  MPI_Comm comm;
  PetscObjectGetComm((PetscObject)p_ref, &comm);
  
  Pot pot; PotCreate(comm, &pot);
  *p_ierr = PotSetFromStr(pot, str);
  
  PetscScalar y_ref[1];
  PetscScalar y[1];
  PetscScalar x[1] = {x0};
  PFApply(p_ref, 1, x, y_ref);
  PFApply(pot,   1, x, y);
  ASSERT_DOUBLE_EQ(PetscRealPart(y_ref[0]), PetscRealPart(y[0]));
  PFDestroy(&pot);

}
TEST(TestPOT, str2pot) {

  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;

  char str_sto[] = "sto 1.1 2 1.2";
  Pot p_sto; PotCreate(comm, &p_sto); PotSetSlater(p_sto, 1.1, 2, 1.2);
  test_str2pot(str_sto, p_sto, 0.2, &ierr);

  char str_rm[] = "pow 1.2 -3";
  Pot p_rm; PotCreate(comm, &p_rm); PotSetPower(p_rm, 1.2, -3);
  test_str2pot(str_rm, p_rm, 1.1, &ierr); 

  char str_comb[] = " sto 1.1 2 1.2 |  pow 1.2 -3 ";
  Pot pfs[2]; pfs[0] = p_sto; pfs[1] = p_rm;
  Pot p_comb; PotCreate(comm, &p_comb); PotSetCombination(p_comb, 2, pfs);
  test_str2pot(str_comb, p_comb, 1.3, &ierr);
  PFDestroy(&p_comb);

  PFDestroy(&p_sto);
  PFDestroy(&p_rm);

  Pot p; PotCreate(comm, &p);
  /* below comments produce error
     PotSetFromOptions2(p, "driv");
  ierr = PotSetFromStr(p, "sto 0 -2 3"); 
  ierr = PotSetFromStr(p, "sto 2.0 -2"); 
  ierr = PotSetFromStr(p, "pow 0.0 2"); 
  ierr = PotSetFromStr(p, "pow 1.1"); 
  */
  PFDestroy(&p);


}
TEST(TestPOT, sum_pot) {
  
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;

  Pot vs[2]; 
  PotCreate(comm, &vs[0]); PotSetSlater(vs[0], 1.1, 2, 0.3);
  PotCreate(comm, &vs[1]); PotSetSlater(vs[1], 0.2, 3, 1.3);

  Pot pot; PotCreate(comm, &pot);
  ierr = PotSetCombination(pot, 2, vs);

  PetscScalar x[1] = {0.55};
  PetscScalar y1[1], y2[1], y_calc[1];
  PFApply(vs[0], 1, x, y1);
  PFApply(vs[1], 1, x, y2);
  PFApply(pot,   1, x, y_calc);
  
  ASSERT_DOUBLE_EQ(PetscRealPart(y1[0]+y2[0]), PetscRealPart(y_calc[0]));

  //  PFView(pot, PETSC_VIEWER_STDOUT_SELF);

}
TEST(TestPOT, product) {
  MPI_Comm comm = MPI_COMM_SELF;
  
  Pot vs[2];
  PotCreate(comm, &vs[0]); PotSetSlater(vs[0], 1.1, 2, 0.2);
  PotCreate(comm, &vs[1]); PotSetPower(vs[1], 1.3, 3);
  
  Pot pot;
  PotCreate(comm, &pot);     
  PotSetProduct(pot, 2, vs);

  PetscScalar x[3] = {0.55, 0.12, 1.1};
  PetscScalar y0[3], y1[3], y[3];
  PFApply(vs[0], 3, x, y0);
  PFApply(vs[1], 3, x, y1);
  PFApply(pot,   3, x, y);

  //PFView(pot, PETSC_VIEWER_STDOUT_SELF);
  
  for(int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(PetscRealPart(y0[i]*y1[i]),
		     PetscRealPart(y[i])) << i;
  
}
TEST(TestPOT, rbessel) {
  MPI_Comm comm = MPI_COMM_SELF;  
  Pot pot; 
  PotCreate(comm, &pot);
  double k = 1.1;
  PotSetRBessel(pot, 0, k);
  PetscScalar x[3] = {0.55, 0.12, 1.1};
  PetscScalar y[3];
  PFApply(pot, 3, x, y);
  for(int i = 0; i < 3; i++)
    EXPECT_DOUBLE_EQ(PetscRealPart(y[i]),
		     PetscRealPart(sin(k*x[i])));

}
TEST(TestPOT, Harmonic) {
  Pot harm; PotCreate(MPI_COMM_SELF, &harm); PotSetHarm(harm, 2.5);

  if(getenv("SHOW_DEBUG"))
    PFView(harm, PETSC_VIEWER_STDOUT_SELF);

  PetscScalar y[1];
  PetscScalar x[1] = {0.2};
  PFApply(harm, 1, x, y);
  ASSERT_DOUBLE_EQ(2.5*0.5*0.2*0.2, PetscRealPart(y[0]));
  PFDestroy(&harm);
}
TEST(TestPOT, Slater) {
  Pot slater; PotCreate(MPI_COMM_SELF, &slater); 
  PotSetSlater(slater, 2.5, 2, 3.1);

  if(getenv("SHOW_DEBUG"))
    PFView(slater, PETSC_VIEWER_STDOUT_SELF);

  PetscScalar y[1];
  PetscScalar x[1] = {0.2};
  PFApply(slater, 1, x, y);
  ASSERT_DOUBLE_EQ(2.5*0.2*0.2*exp(-3.1*0.2), PetscRealPart(y[0]));
  PFDestroy(&slater);
}

int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}
