#include <slepceps.h>
#include <gtest/gtest.h>
#include <rescol/cscaling.h>

static char help[] = "Unit test for scale.c";

TEST(TestScaling, Uniform) {

  MPI_Comm comm = MPI_COMM_SELF;
  PetscReal theta = M_PI*10.0/18.0;
  CScaling cscaling; CScalingCreate(comm, &cscaling);
  CScalingSetUniformCS(cscaling, theta);

  if(getenv("SHOW_DEBUG"))
    CScalingView(cscaling, PETSC_VIEWER_STDOUT_SELF);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 0.2; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;
  ierr = CScalingCalc(cscaling, xs, 10, qr, Rr); ASSERT_EQ(0, ierr);  

  double eps = 0.00000000001;
#if defined(PETSC_USE_COMPLEX)
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 0.4*cos(theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(Rr[2]), 0.4*sin(theta), eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), cos(theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(qr[2]), sin(theta), eps);
#else
  ASSERT_NEAR(Rr[0], 0.0, eps);
  ASSERT_NEAR(Rr[2], 0.4*cos(theta), eps);
  ASSERT_NEAR(qr[0], cos(-theta), eps);
  ASSERT_NEAR(qr[2], cos(-theta), eps);
#endif

  ierr = CScalingDestroy(&cscaling); ASSERT_EQ(0, ierr);
}
TEST(TestScaling, SharpECS) {
  
  MPI_Comm comm = MPI_COMM_SELF;
  PetscReal theta = M_PI*10.0/180.0;
  PetscReal r0 = 4.0;
  CScaling cscaling; CScalingCreate(comm, &cscaling);
  CScalingSetSharpECS(cscaling, r0, theta);

  if(getenv("SHOW_DEBUG"))
    CScalingView(cscaling, PETSC_VIEWER_STDOUT_SELF);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 1.0; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;
  ierr = CScalingCalc(cscaling, xs, 10, qr, Rr);  ASSERT_EQ(0, ierr);  

  double eps = 0.00000000001;

#if defined(PETSC_USE_COMPLEX)
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 2.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[5]), 4.0 + 1.0*cos(theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(Rr[5]), 1.0*sin(theta), eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[2]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[6]), cos(theta), eps);
  ASSERT_NEAR(PetscImaginaryPart(qr[6]), sin(theta), eps);
#else
  ASSERT_NEAR(Rr[0], 0.0, eps);
  ASSERT_NEAR(Rr[2], 2.0, eps);
  ASSERT_NEAR(Rr[5], 4.0 + 1.0*cos(theta), eps);
  ASSERT_NEAR(qr[0], 1.0, eps);
  ASSERT_NEAR(qr[2], 1.0, eps);
  ASSERT_NEAR(qr[6], cos(theta), eps);
#endif  

  ierr = CScalingDestroy(&cscaling); ASSERT_EQ(0, ierr);  
}
TEST(TestScaling, none) {
  
  MPI_Comm comm = MPI_COMM_SELF;
  CScaling cscaling; CScalingCreate(comm, &cscaling); CScalingSetNone(cscaling);

  if(getenv("SHOW_DEBUG"))
    CScalingView(cscaling, PETSC_VIEWER_STDOUT_SELF);

  PetscReal xs[10];
  PetscScalar Rr[10];
  PetscScalar qr[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 1.1; Rr[i] = 0.0; qr[i] = 0.0;
  }

  PetscErrorCode ierr;
  ierr = CScalingCalc(cscaling, xs, 10, qr, Rr);  ASSERT_EQ(0, ierr);  

  double eps = 0.00000000001;
 
  ASSERT_NEAR(PetscRealPart(Rr[0]), 0.0, eps);
  ASSERT_NEAR(PetscRealPart(Rr[2]), 2.2, eps);
  ASSERT_NEAR(PetscRealPart(Rr[5]), 5.5, eps);
  ASSERT_NEAR(PetscRealPart(qr[0]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[2]), 1.0, eps);
  ASSERT_NEAR(PetscRealPart(qr[6]), 1.0, eps);

  ierr = CScalingDestroy(&cscaling); ASSERT_EQ(0, ierr);  

}
TEST(TestScaling , Copy) {
  
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  PetscReal theta = M_PI*10.0/180.0;
  PetscReal r0 = 4.0;
  CScaling cscaling;
  ierr = CScalingCreate(comm, &cscaling); ASSERT_EQ(0, ierr);
  ierr = CScalingSetSharpECS(cscaling, r0, theta); ASSERT_EQ(0, ierr);

  CScaling cscaling2;
  ierr = CScalingCreate(comm, &cscaling2); ASSERT_EQ(0, ierr);
  ierr = CScalingCopy(cscaling, cscaling2); ASSERT_EQ(0, ierr);

  if(getenv("SHOW_DEBUG"))
    CScalingView(cscaling, PETSC_VIEWER_STDOUT_SELF);

  PetscReal xs[10];
  PetscScalar Rr[10], Rr2[10];
  PetscScalar qr[10], qr2[10];
  for(int i = 0; i < 10; i++) {
    xs[i] = i * 1.0; Rr[i] = 0.0; qr[i] = 0.0;
  }

  ierr = CScalingCalc(cscaling, xs, 10, qr, Rr);  ASSERT_EQ(0, ierr);  
  ierr = CScalingCalc(cscaling2,xs, 10, qr2, Rr2);  ASSERT_EQ(0, ierr);  

  double eps = 0.00000000001;

  for(int i = 0; i < 10; i++) {
    ASSERT_NEAR(PetscRealPart(Rr[i]), PetscRealPart(Rr2[i]), eps);
    ASSERT_NEAR(PetscRealPart(qr[i]), PetscRealPart(qr2[i]), eps);
    ASSERT_NEAR(PetscImaginaryPart(Rr[i]), PetscImaginaryPart(Rr2[i]), eps);
    ASSERT_NEAR(PetscImaginaryPart(qr[i]), PetscImaginaryPart(qr2[i]), eps);
  }

  CScalingDestroy(&cscaling);
  CScalingDestroy(&cscaling2);
  
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

