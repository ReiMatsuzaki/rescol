#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/scale.h>

static char help[] = "Unit test for scale.c";

TEST(TestScaling, Uniform) {

  MPI_Comm comm = MPI_COMM_SELF;
  PetscReal theta = M_PI*10.0/18.0;
  Scaler scaler;
  ScalerCreateUniform(&scaler, comm, theta);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 0.2; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;
  ierr = ScalerSetRr(scaler, xs, 10, Rr); ASSERT_EQ(0, ierr);  
  ierr = ScalerSetQr(scaler, xs, 10, qr); ASSERT_EQ(0, ierr);

  double eps = 0.00000000001;
#if defined(PETSC_USE_COMPLEX)
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 0.4*cos(-theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(Rr[2]), 0.4*sin(-theta), eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), cos(-theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(qr[2]), sin(-theta), eps);
#else
  ASSERT_NEAR(Rr[0], 0.0, eps);
  ASSERT_NEAR(Rr[2], 0.4*cos(-theta), eps);
  ASSERT_NEAR(qr[0], cos(-theta), eps);
  ASSERT_NEAR(qr[2], cos(-theta), eps);
#endif
  
  

  ierr = ScalerDestroy(&scaler); ASSERT_EQ(0, ierr);
}
TEST(TestScaling, SharpECS) {
  
  MPI_Comm comm = MPI_COMM_SELF;
  PetscReal theta = M_PI*10.0/180.0;
  PetscReal r0 = 4.0;
  Scaler scaler;
  ScalerCreateSharpECS(&scaler, comm, r0, theta);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 1.0; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;

  ierr = ScalerSetRr(scaler, xs, 10, Rr); ASSERT_EQ(0, ierr);  
  ierr = ScalerSetQr(scaler, xs, 10, qr); ASSERT_EQ(0, ierr);

  double eps = 0.00000000001;

#if defined(PETSC_USE_COMPLEX)
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 2.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[5]), 4.0 + 1.0*cos(-theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(Rr[5]), 1.0*sin(-theta), eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[2]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[6]), cos(-theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(qr[6]), sin(-theta), eps);
#else
  ASSERT_NEAR(Rr[0], 0.0, eps);
  ASSERT_NEAR(Rr[2], 2.0, eps);
  ASSERT_NEAR(Rr[5], 4.0 + 1.0*cos(-theta), eps);
  ASSERT_NEAR(qr[0], 1.0, eps);
  ASSERT_NEAR(qr[2], 1.0, eps);
  ASSERT_NEAR(qr[6], cos(-theta), eps);
#endif  

  ierr = ScalerDestroy(&scaler); ASSERT_EQ(0, ierr);  
}
TEST(TestScaling, none) {
  
  MPI_Comm comm = MPI_COMM_SELF;
  Scaler scaler;
  ScalerCreateNone(&scaler, comm);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 1.1; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;

  ierr = ScalerSetRr(scaler, xs, 10, Rr); ASSERT_EQ(0, ierr);  
  ierr = ScalerSetQr(scaler, xs, 10, qr); ASSERT_EQ(0, ierr);

  double eps = 0.00000000001;
  
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 2.2, eps);
  ASSERT_NEAR(PetscRealPart(Rr[5]), 5.5, eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[2]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[6]), 1.0, eps);

  ierr = ScalerDestroy(&scaler); ASSERT_EQ(0, ierr);  

}
int main (int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help);
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
  SlepcFinalize();
  return 0;
}

