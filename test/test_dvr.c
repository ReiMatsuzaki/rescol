#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/mat.h>
#include <rescol/dvr.h>

int dgetrf_(long*, long*, double*, long*, long*, long*);
int dgetri_(long*, double*, long*, long*, double*, long*, long*);

static char help[] = "Unit test for dvr.c \n\n";

int testXS() {
  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 3;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  ASSERT_DOUBLE_EQ(0.0, dvr->xs[0]);
  ASSERT_DOUBLE_EQ(0.5, dvr->xs[1]);  
  ASSERT_DOUBLE_EQ(1.0, dvr->xs[2]);

  ASSERT_DOUBLE_NEAR(0.1666666666, dvr->ws[0], 0.000000001);
  ASSERT_DOUBLE_NEAR(0.6666666666, dvr->ws[1], 0.000000001);
  ASSERT_DOUBLE_NEAR(0.1666666666, dvr->ws[2], 0.000000001);

  ASSERT_DOUBLE_EQ(0.5, dvr->xs_basis[0]);
  ASSERT_DOUBLE_EQ(1.0, dvr->xs_basis[1]);
  ASSERT_DOUBLE_EQ(1.5, dvr->xs_basis[2]);
  ASSERT_DOUBLE_EQ(2.0, dvr->xs_basis[3]);
  ASSERT_DOUBLE_EQ(2.5, dvr->xs_basis[4]);  

  DVRDestroy(&dvr);
  return 0;
}
int testSR1LSMat() {
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  
  int nq = 3;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat S; DVRCreateR1LSMat(dvr, &S); DVRSR1LSMat(dvr, S);
  
  PetscInt n, m;
  MatGetSize(S, &n, &m);
  ASSERT_EQ(15, n); ASSERT_EQ(15, m);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  ierr = MatGetRow(S, 0, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_DOUBLE_EQ(dvr->ws[0], row[0]);
  ASSERT_DOUBLE_EQ(0.1+0.2/3.0, row[0]);
  ASSERT_EQ(0, cols[0]);
  ierr = MatRestoreRow(S, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr = MatGetRow(S, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_DOUBLE_EQ(dvr->ws[1], row[0]);
  ASSERT_EQ(1, cols[0]);
  ierr = MatRestoreRow(S, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  DVRDestroy(&dvr);
  MatDestroy(&S);   
  return 0;

}
int testD2R1LSMat() {
  /*
    1: 2.66666667 -5.33333333  2.66666667
   */
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 3;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat D; DVRCreateR1LSMat(dvr, &D); DVRD2R1LSMat(dvr, D);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;

  ierr = MatGetRow(D, 0, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(3, ncols);
  ASSERT_EQ(0, cols[0]); ASSERT_EQ(1, cols[1]); ASSERT_EQ(2, cols[2]);
  ASSERT_DOUBLE_NEAR(-2.3333333333333, row[0], 0.0000000001); 
  ASSERT_DOUBLE_NEAR(2.66666666666666, row[1], 0.0000000001);
  ASSERT_DOUBLE_NEAR(-0.33333333333333, row[2], 0.0000000001);
  ierr = MatRestoreRow(D, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr = MatGetRow(D, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(3, ncols);
  ASSERT_EQ(0, cols[0]); ASSERT_EQ(1, cols[1]); ASSERT_EQ(2, cols[2]);
  ASSERT_DOUBLE_EQ(8.0/3.0, row[0]); 
  ASSERT_DOUBLE_EQ(-16.0/3.0, row[1]);
  ASSERT_DOUBLE_EQ(8.0/3.0, row[2]);
  ierr = MatRestoreRow(D, 1, &ncols, &cols, &row); CHKERRQ(ierr);

  DVRDestroy(&dvr);
  MatDestroy(&D);
  return 0;
}
int testENR1LSMat() {

  /*
    |...|...|...|... |
   */
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 3;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat S; DVRCreateR1LSMat(dvr, &S); DVRSR1LSMat(dvr, S);


  Mat V; DVRCreateR1LSMat(dvr, &V); DVRENR1LSMat(dvr, 0, 0.0, V);
  
  PetscInt n, m;
  MatGetSize(V, &n, &m);
  ASSERT_EQ(15, n); ASSERT_EQ(15, m);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  ierr = MatGetRow(V, 0, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_DOUBLE_EQ(0.0, row[0]);
  ASSERT_EQ(0, cols[0]);
  ierr = MatRestoreRow(V, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr = MatGetRow(V, 1, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(1, ncols);
  ASSERT_DOUBLE_EQ(dvr->ws[1]/dvr->xs[1], row[0]);
  ASSERT_EQ(1, cols[0]);
  ierr = MatRestoreRow(V, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  DVRDestroy(&dvr);
  MatDestroy(&S);
  MatDestroy(&V);
  return 0;
}
int testENR1Mat() {

  /*
    |...|...|...|... |
   */

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 3;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 5);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);
  
  Mat V; DVRCreateR1Mat(dvr, &V); DVRENR1Mat(dvr, 0, 0.0, V);
  Mat D; DVRCreateR1LSMat(dvr, &D); DVRD2R1LSMat(dvr, D);
  Mat S; DVRCreateR1LSMat(dvr, &S); DVRSR1LSMat(dvr, S);
  
  PetscInt n, m;
  MatGetSize(V, &n, &m);

  ASSERT_EQ(7, n); ASSERT_EQ(7, m);

  MatDestroy(&V);
  MatDestroy(&D);
  MatDestroy(&S);
  DVRDestroy(&dvr);

  return 0;
}
int testLSMat_to_Mat() {
  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 4;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 10.0, 4);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat S; DVRCreateR1LSMat(dvr, &S);
  int ne; BPSGetNumEle(bps, &ne);
  for(int i = 0; i < 12; i++)
    for(int j = 0; j < 12; j++) {
      PetscScalar v = i/nq+ne*(j/nq);
      MatSetValue(S, j, i, v, INSERT_VALUES);
    }
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);

  Mat SS;
  DVRLSMatToMat(dvr, S, MAT_INITIAL_MATRIX, &SS);

  PetscInt n, m;
  MatGetSize(SS, &n, &m);
  ASSERT_EQ(8, n); ASSERT_EQ(8, m);

  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;

  ierr = MatGetRow(SS, 0, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_DOUBLE_NEAR(0.0, row[0], 0.0000000001); 
  ierr = MatRestoreRow(SS, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr = MatGetRow(SS, 5, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_DOUBLE_NEAR(20.0, row[2], 0.0000000001);
  ierr = MatRestoreRow(SS, 5, &ncols, &cols, &row); CHKERRQ(ierr);

  return 0;
}
int testHAtom() {

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 5;
  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 20.0, 20, 5.0);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat H; DVRCreateR1LSMat(dvr, &H); DVRD2R1LSMat(dvr, H); MatScale(H, -0.5);
  Mat V; DVRCreateR1LSMat(dvr, &V); DVRENR1LSMat(dvr, 0, 0.0, V);
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  Mat S; DVRCreateR1LSMat(dvr, &S); DVRSR1LSMat(dvr, S);

  Mat HH; DVRLSMatToMat(dvr, H, MAT_INITIAL_MATRIX, &HH);
  Mat SS; DVRLSMatToMat(dvr, S, MAT_INITIAL_MATRIX, &SS);

  EPS eps; 
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, HH, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSJD);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6);
  EPSSolve(eps);

  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  DVRDestroy(&dvr);
  MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  MatDestroy(&HH); MatDestroy(&SS);
  EPSDestroy(&eps);

  return 0;
}
int testHAtom2() {

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 5;
  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 20.0, 20, 5.0);
  DVR dvr; DVRCreate(comm, &dvr); DVRSetKnots(dvr, nq, bps); DVRSetUp(dvr);

  Mat H; DVRCreateR1Mat(dvr, &H); DVRD2R1Mat(dvr, H); MatScale(H, -0.5);
  Mat V; DVRCreateR1Mat(dvr, &V); DVRENR1Mat(dvr, 0, 0.0, V);
  Mat L; DVRCreateR1Mat(dvr, &L); DVRR2invR1Mat(dvr, L);

  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(H, 1*(1+1)/2.0, L, DIFFERENT_NONZERO_PATTERN);

  EPS eps; 
  EPSCreate(PETSC_COMM_SELF, &eps);
  EPSSetOperators(eps, H, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSJD);
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6);
  EPSSolve(eps);

  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-1.0/8.0, kr, pow(10.0, -5.0));

  DVRDestroy(&dvr);
  MatDestroy(&H);
  MatDestroy(&V);
  MatDestroy(&L);
  EPSDestroy(&eps);

  return 0;
}
int testHeAtom() {
  ASSERT_EQ(2, 1);
/*

  PetscErrorCode ierr;

  PrintTimeStamp(PETSC_COMM_SELF, "Init", NULL);
  DVR dvr;
  int nq = 5;
  double rmax = 20.0;
  int num_zs = 21; // number of R2 basis is (5*20)**2=10000
  BPS bps; BPSCreate(&bps, PETSC_COMM_SELF); BPSSetExp(bps, rmax, num_zs, 5.0);
  DVRCreate(&dvr, nq, bps, PETSC_COMM_SELF);

  Mat S_R1, T_R1, V_R1;  
  PrintTimeStamp(dvr->comm, "R1", NULL);
  ierr = DVRSetSR1Mat(dvr, &S_R1); CHKERRQ(ierr);
  ierr = DVRSetD2R1Mat(dvr, &T_R1); CHKERRQ(ierr); MatScale(T_R1, -0.5);  
  ierr = DVRSetENR1Mat(dvr, 0, 0.0, &V_R1); CHKERRQ(ierr); MatScale(V_R1, -1.0); 

  Mat H;
  PrintTimeStamp(dvr->comm, "T", NULL);
  ierr = MatSetSynthesizeFast(S_R1, T_R1, dvr->comm, &H); CHKERRQ(ierr);

  Mat T2;
  ierr = MatSetSynthesizeFast(T_R1, S_R1, dvr->comm, &T2); CHKERRQ(ierr);
  ierr = MatAXPY(H, 1.0, T2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&T2);

  Mat V1;
  PrintTimeStamp(dvr->comm, "e-n", NULL);
  ierr = MatSetSynthesizeFast(S_R1, V_R1, dvr->comm, &V1); CHKERRQ(ierr);
  ierr = MatAXPY(H, 1.0, V1, SUBSET_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V1);

  Mat V2;
  ierr = MatSetSynthesizeFast(V_R1, S_R1, dvr->comm, &V2); CHKERRQ(ierr);
  ierr = MatAXPY(H, 1.0, V2, SUBSET_NONZERO_PATTERN); CHKERRQ(ierr);
  MatDestroy(&V2);

  Mat Vee;
  PrintTimeStamp(dvr->comm, "e-e", NULL);
  ierr = DVRSetEER2Mat(dvr, 0, &Vee); CHKERRQ(ierr);
  ierr = MatAXPY(H, 1.0, Vee, SUBSET_NONZERO_PATTERN); CHKERRQ(ierr); 
  MatDestroy(&Vee); 

  EPS eps; 
  PrintTimeStamp(dvr->comm, "eps", NULL);
  ierr = EPSCreate(PETSC_COMM_SELF, &eps); CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, H, NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
  ierr = EPSSetTarget(eps, -3.0); CHKERRQ(ierr);
  PrintTimeStamp(dvr->comm, "solve", NULL);
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  int nconv;
  PetscScalar kr;
  PrintTimeStamp(dvr->comm, "output", NULL);
  ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
  ASSERT_TRUE(nconv > 0);
  ierr = EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL); CHKERRQ(ierr);
  ASSERT_DOUBLE_NEAR(-1.0/8.0, kr, pow(10.0, -5.0));
  */

  return 0;  
}
int testLapack() {
  long m = 2;
  long n = m;
  long lda = m;
  double A[2*2];
  A[0] = 2.0; A[2] = 3.0;
  A[1] = 1.0; A[3] = 0.5;
  long info;
  long ipiv[2];
  long lwork = 2;
  double work[2];

  dgetrf_(&m, &n, A, &lda, ipiv, &info);
  dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);

  ASSERT_DOUBLE_EQ(-0.25, A[0]);
  ASSERT_DOUBLE_EQ(0.5, A[1]);  
  ASSERT_DOUBLE_EQ(1.5, A[2]);  
  ASSERT_DOUBLE_EQ(-1.0, A[3]);    


  return 0;
}

int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);
  PrintTimeStamp(PETSC_COMM_SELF, "testXS", NULL);
  testXS();
  PrintTimeStamp(PETSC_COMM_SELF, "testSR1LSMat", NULL);
  testSR1LSMat();
  PrintTimeStamp(PETSC_COMM_SELF, "testD2R1LSMat", NULL);
  testD2R1LSMat();
  PrintTimeStamp(PETSC_COMM_SELF, "testENR1LSMat", NULL);
  testENR1LSMat();
  testENR1Mat();
  
  // testLSMat_to_Mat();
  PrintTimeStamp(PETSC_COMM_SELF, "testHAtom", NULL);
  testHAtom();
  PrintTimeStamp(PETSC_COMM_SELF, "testHAtom2", NULL);
  testHAtom2();


  //testHeAtom();
  //testLapack();
  SlepcFinalize();
  return 0;
}

