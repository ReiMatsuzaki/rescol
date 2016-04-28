#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/mat.h>
#include <rescol/dvr.h>


int dgetrf_(long*, long*, double*, long*, long*, long*);
int dgetri_(long*, double*, long*, long*, double*, long*, long*);

static char help[] = "Unit test for dvr.c \n\n";

int testLS() {

  PetscScalar xs_c[6] = {1.1, 1.2, 1.4, 1.4, 1.9, 2.1};
  PetscScalar xs_r[6] = {1.1, 1.2, 1.4, 1.4, 1.9, 2.1};
  PetscReal zs[3]   = {1.1,        1.4,           2.1};
  
  PetscScalar y;
  for(int iele = 0; iele < 2; iele++) {
    for(int iq = 0; iq < 3; iq++) {
      for(int ils = 0; ils < 3; ils++) {
	ValueLS(xs_c, zs, 2, 3, iele, ils,
		xs_r[iele*3+iq],
		xs_c[iele*3+iq], &y);
	if(ils == iq) {
	  ASSERT_SCALAR_EQ(1.0, y);
	}
	else {
	  ASSERT_SCALAR_EQ(0.0, y);
	}
      }
    }
  }
  return 0;
}
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
int testECSMatrix() {

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 5;
  BPS bps;
  BPSCreate(comm, &bps);
  BPSSetLine(bps, 10.0, 11);
  DVR dvr;
  DVRCreate(comm, &dvr);
  DVRSetKnots(dvr, nq, bps);
  DVRSetCScaling(dvr, 5.0, 20.0);
  
  DVRSetUp(dvr);

  //  DVRView(dvr, PETSC_VIEWER_STDOUT_SELF);

  Mat S;
  DVRCreateR1LSMat(dvr, &S);
  DVRSR1LSMat(dvr, S);
  Mat T;
  DVRCreateR1LSMat(dvr, &T);
  DVRD2R1LSMat(dvr, T);
  MatScale(T, -0.5);
  Mat V;
  DVRCreateR1LSMat(dvr, &V);
  DVRENR1LSMat(dvr, 2, 0.0, V);
  Mat V0;
  DVRCreateR1Mat(dvr, &V0);
  DVRENR1Mat(dvr, 2, 0.0, V0);

  Mat SS;
  DVRLSMatToMat(dvr, S, MAT_INITIAL_MATRIX, &SS);  
  Mat TT;
  DVRLSMatToMat(dvr, T, MAT_INITIAL_MATRIX, &TT);
  Mat VV;
  DVRLSMatToMat(dvr, V, MAT_INITIAL_MATRIX, &VV);  

  // MatView(V, PETSC_VIEWER_STDOUT_SELF);
  // MatView(dvr->T, PETSC_VIEWER_STDOUT_SELF);
  // MatView(dvr->TT, PETSC_VIEWER_STDOUT_SELF);
  // MatView(VV,PETSC_VIEWER_STDOUT_SELF);

  //  ASSERT_MAT_EQ(V0, VV);
  test_mat_eq(V0, VV, pow(10.0, -10.0), __FILE__, __LINE__);

  int n; DVRGetSize(dvr, &n);
  PetscScalar *s;  PetscMalloc1(n*n, &s);
  PetscScalar *t;  PetscMalloc1(n*n, &t);
  PetscScalar *v;  PetscMalloc1(n*n, &v);
  int *idx;        PetscMalloc1(n,   &idx);
  for(int i = 0; i < n; i++) { idx[i] = i; }
  
  MatGetValues(SS, n, idx, n, idx, s);
  MatGetValues(TT, n, idx, n, idx, t);
  MatGetValues(VV, n, idx, n, idx, v);

  for(int i = 0; i < n; i++) {
    PetscScalar xi = dvr->xs_basis_c[i];
    //printf("%f, %f\n", creal(xi), cimag(xi));
    for(int j = 0; j < n; j++) {
      ASSERT_SCALAR_EQ(i==j?1.0:0.0,    s[n*i+j]);
      ASSERT_SCALAR_EQ(i==j?1.0/xi:0.0, v[n*i+j]);
    }
  }

  PetscFree(s); PetscFree(t); PetscFree(v);
  PetscFree(idx);
  MatDestroy(&S); MatDestroy(&T); MatDestroy(&V); MatDestroy(&V0);
  MatDestroy(&SS);MatDestroy(&TT);MatDestroy(&VV);
  DVRDestroy(&dvr);
  return 0;
}
int testECSTMat() {

  PetscErrorCode ierr;

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 5;
  BPS bps;
  ierr = BPSCreate(comm, &bps); CHKERRQ(ierr);
  ierr = BPSSetLine(bps, 10.0, 11); CHKERRQ(ierr);

  DVR dvr_c;
  ierr = DVRCreate(comm, &dvr_c); CHKERRQ(ierr);
  ierr = DVRSetKnots(dvr_c, nq, bps); CHKERRQ(ierr);
  ierr = DVRSetCScaling(dvr_c, 0.0, 20.0); CHKERRQ(ierr);
  ierr = DVRSetUp(dvr_c); CHKERRQ(ierr);

  BPS bps2;
  ierr = BPSCreate(comm, &bps2); CHKERRQ(ierr);
  ierr = BPSSetLine(bps2, 10.0, 11); CHKERRQ(ierr);

  DVR dvr_r;
  ierr = DVRCreate(comm, &dvr_r); CHKERRQ(ierr);
  ierr = DVRSetKnots(dvr_r, nq, bps2); CHKERRQ(ierr);
  ierr = DVRSetUp(dvr_r); CHKERRQ(ierr);
  
  Mat T_r;
  ierr = DVRCreateR1Mat(dvr_r, &T_r); CHKERRQ(ierr);
  ierr = DVRD2R1Mat(dvr_r, T_r); CHKERRQ(ierr);
  
  Mat T_c;
  ierr = DVRCreateR1Mat(dvr_c, &T_c); CHKERRQ(ierr);
  ierr = DVRD2R1Mat(dvr_c, T_c); CHKERRQ(ierr);

  Mat V_r;
  ierr = DVRCreateR1Mat(dvr_r, &V_r); CHKERRQ(ierr);
  ierr = DVRENR1Mat(dvr_r, 0, 0.0, V_r); CHKERRQ(ierr);

  Mat V_c;
  ierr = DVRCreateR1Mat(dvr_c, &V_c); CHKERRQ(ierr);
  ierr = DVRENR1Mat(dvr_c, 0, 0.0, V_c); CHKERRQ(ierr);  

  PetscScalar scale_factor = cexp(-I*20.0*M_PI/180.0);
  ierr = MatScale(T_r, scale_factor*scale_factor); CHKERRQ(ierr);
  ierr = MatScale(V_r, scale_factor); CHKERRQ(ierr);

  test_mat_eq(T_r, T_c, pow(10.0, -10.0), __FILE__, __LINE__);
  test_mat_eq(V_r, V_c, pow(10.0, -10.0), __FILE__, __LINE__);

  DVRDestroy(&dvr_r);
  DVRDestroy(&dvr_c);
  MatDestroy(&T_r);
  MatDestroy(&T_c);
  MatDestroy(&V_r);
  MatDestroy(&V_c);

  return 0;

}
int testECSVector() {

  MPI_Comm comm = PETSC_COMM_SELF;
  int nq = 3;
  BPS bps;
  BPSCreate(comm, &bps);
  BPSSetLine(bps, 2.0, 3);
  DVR dvr;
  DVRCreate(comm, &dvr);
  DVRSetKnots(dvr, nq, bps);
  DVRSetCScaling(dvr, 1.0, 20.0);
  
  DVRSetUp(dvr);
  //  printf("basis size: %d\n", dvr->num_basis);
  //  DVRView(dvr, PETSC_VIEWER_STDOUT_SELF);

  Pot slater;
  PotCreate(comm, &slater);
  PotSetSlater(slater, 1.1, 2, 1.2);

  Vec m_ls;
  DVRCreateR1LSVec(dvr, &m_ls);
  DVRPotR1LSVec(dvr, slater, m_ls);

  Vec m_trans;
  DVRCreateR1Vec(dvr, &m_trans);
  DVRLSVecToVec(dvr, m_ls, m_trans);
  
  Vec m_basis;
  DVRCreateR1Vec(dvr, &m_basis);
  DVRPotR1Vec(dvr, slater, m_basis);

  /*
  VecView(m_ls, PETSC_VIEWER_STDOUT_SELF);
  VecView(m_trans, PETSC_VIEWER_STDOUT_SELF);
  VecView(m_basis, PETSC_VIEWER_STDOUT_SELF);
  MatView(dvr->TT, PETSC_VIEWER_STDOUT_SELF);
  */

/*
  PetscScalar *xs_trans, *xs_basis;
  VecGetArray(m_trans, &xs_trans);
  VecGetArray(m_basis, &xs_basis);

  int n; DVRGetSize(dvr, &n);
  for(int i = 0; i < n; i++) {
    ASSERT_SCALAR_EQ(xs_trans[i], xs_basis[i]);
  }

  VecRestoreArray(m_trans, &xs_trans);
  VecRestoreArray(m_basis, &xs_basis);
*/

  test_vec_eq(m_trans, m_basis, __FILE__, __LINE__);

  VecDestroy(&m_ls); VecDestroy(&m_trans); VecDestroy(&m_basis);
  DVRDestroy(&dvr);
  PFDestroy(&slater);
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

  PrintTimeStamp(PETSC_COMM_SELF, "test_valueLS", NULL);
  testLS();

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

  PrintTimeStamp(PETSC_COMM_SELF, "testECSMatrix", NULL);
  testECSMatrix();
  PrintTimeStamp(PETSC_COMM_SELF, "testECS_T", NULL);
  testECSTMat();
  PrintTimeStamp(PETSC_COMM_SELF, "testECSVector", NULL);
  testECSVector();

  //testHeAtom();
  //testLapack();
  SlepcFinalize();
  return 0;
}

