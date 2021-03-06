#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include "../include/math.h"
#include "../include/bspline.h"
#include "../include/pot.h"
#include "../include/eeps.h"

static char help[] = "Unit test for bspline.c \n\n";

int testNumBSpline() {
  PrintTimeStamp(PETSC_COMM_SELF, "num bs", NULL);
  /*
    order = 4
----|. . . . . |   0 (not used)
 ---|- . . . . |   1 (not used)
  --|- - . . . |   2 
   -|- - - . .|    3 
    |- - - - .|    4
    |. - - - -|    5
    |. . - - -|-    6
    |. . . - -|--    7    
    |. . . . -|---    8   
    
    order = 3
    num_ele = 5 case
 ---|. . . . . |   0 (not used)
  --|- . . . . |   1 (not used)
   -|- - . . .|    2 
    |- - - . .|    3
    |. - - - .|    4
    |. . - - -|    5
    |. . . - -|-    6

    order = 2
    num_ele = 5 case
 --|. . . . .|    0 (not used)
  -|- . . . .|    1 (not used)
   |- - . . .|    2 (0)
   |. - - . .|    3 (1)
   |. . - - .|    4 (2)
   |. . . - -|    5 (3)
   |. . . . -|-   6 (not used)
   |. . . . .|--  7 (not used)

    where (-) means non zero (.) means zero
   */
  ASSERT_EQ(5, NumBSpline(3, 5));
  ASSERT_EQ(4, NumBSpline(2, 5));
  return 0;
}
int testCalcBSpline() {
  PrintTimeStamp(PETSC_COMM_SELF, "calc", NULL);

  // see paper
  int order = 3;
  
  // overlaped knot points list
  double ts_r[10] = {0.0, 0.0, 0.0, 
		   1.0, 2.0, 3.0, 4.0,
		   5.0, 5.0, 5.0};
  PetscScalar ts_s[10];
  for(int i = 0; i < 10; i++)
    ts_s[i] = ts_r[i];
  
  // non overlaped points list
  //double zs[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0}; 

  PetscBool zeroq;
  PetscScalar y = 777.0;
  CalcBSpline(order, ts_r, ts_s, 0, 0.0, 0.0, &y, &zeroq); 
  ASSERT_DOUBLE_EQ(1.0, y); ASSERT_FALSE(zeroq);
  CalcBSpline(order, ts_r, ts_s, 0, 1.0, 1.0, &y, &zeroq);
  ASSERT_DOUBLE_EQ(0.0, y); ASSERT_TRUE(zeroq);
  CalcBSpline(order, ts_r, ts_s, 6, 4.0, 4.0, &y, &zeroq); 
  ASSERT_DOUBLE_EQ(0.0, y); ASSERT_TRUE(zeroq);
  
  PetscScalar x = 0.34;
  CalcBSpline(     order, ts_r, ts_s, 2, x, x, &y, &zeroq);
  ASSERT_DOUBLE_EQ(0.5*x*x, y); ASSERT_FALSE(zeroq);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, x, &y, &zeroq);
  ASSERT_DOUBLE_EQ(x, y); ASSERT_FALSE(zeroq);

  x = 2.44;
  CalcBSpline(     order, ts_r, ts_s, 2, x, x, &y, &zeroq);
  ASSERT_DOUBLE_EQ(0.5*x*x-3*x+4.5, y); ASSERT_FALSE(zeroq);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, x, &y, &zeroq);
  ASSERT_DOUBLE_EQ(x-3.0, y); ASSERT_FALSE(zeroq);

  x = 3.44;
  CalcBSpline(     order, ts_r, ts_s, 2, x, x, &y, &zeroq);
  ASSERT_DOUBLE_EQ(0.0, y); ASSERT_TRUE(zeroq);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, x, &y, &zeroq); 
  ASSERT_DOUBLE_EQ(0.0, y); ASSERT_TRUE(zeroq);

  return 0;
}
int testCalcBSplineECS() {
  PrintTimeStamp(PETSC_COMM_SELF, "calc", NULL);

/*
  // overlaped knot points list
  double ts_r[10] = {0.0, 0.0, 0.0, 
		   1.0, 2.0, 3.0, 4.0,
		   5.0, 5.0, 5.0};

  PetscScalar ts_s[10];
  for(int i = 0; i < 10; i++) {
    if(ts_r[i] > 2.5) {
      PetscScalar arg = 0.1I;
      ts_s[i] = 2.5 + (ts_r[i]-2.5)*exp(arg);
    } else {
      ts_s[i] = ts_r[i];
    }
  }
  
  PetscScalar x = 3.1;
  PetscReal h = 0.0001;
  PetscScalar yp, ym, dy;
  
  int order = 8;
  for(int index = 0; index < 8; index++) {
    printf("%d, %f, %f\n", index, creal(yp), cimag(yp));
    CalcBSpline(     order, ts_r, ts_s, index, x+h, x+h, &yp); 
    CalcBSpline(     order, ts_r, ts_s, index, x-h, x-h, &ym); 
    CalcDerivBSpline(order, ts_r, ts_s, index, x, x, &dy); 
    
    ASSERT_DOUBLE_NEAR(creal((yp-ym)/(2.0*h)), creal(dy), 0.001);
    ASSERT_DOUBLE_NEAR(cimag((yp-ym)/(2.0*h)), cimag(dy), 0.001);
  }
  
*/
  return 0;
}
int testNon0QuadIndex() {
  PrintTimeStamp(PETSC_COMM_SELF, "non0", NULL);

  int i0, i1;
  Non0QuadIndex(0, 0, 3, 3*5, &i0, &i1);
  ASSERT_EQ(0, i0); ASSERT_EQ(6, i1);

  Non0QuadIndex(0, 1, 3, 3*5, &i0, &i1);
  ASSERT_EQ(0, i0); ASSERT_EQ(6, i1);

  Non0QuadIndex(1, 0, 3, 3*5, &i0, &i1);
  ASSERT_EQ(0, i0); ASSERT_EQ(6, i1);

  Non0QuadIndex(2, 0, 3, 3*5, &i0, &i1);
  ASSERT_EQ(3, i0); ASSERT_EQ(6, i1);

  Non0QuadIndex(4, 2, 3, 3*5, &i0, &i1);
  ASSERT_EQ(3*3, i0); ASSERT_EQ(3*4, i1);

  Non0QuadIndex(3, 4, 3, 3*5, &i0, &i1);
  ASSERT_EQ(3*3, i0); ASSERT_EQ(3*5, i1);

  Non0QuadIndex(4, 4, 3, 3*5, &i0, &i1);
  ASSERT_EQ(3*3, i0); ASSERT_EQ(3*5, i1);
    
  return 0;
}
int testNon0QuadIndex2() {
  PrintTimeStamp(PETSC_COMM_SELF, "non02", NULL);


  int i0, i1;
  int k  = 2;
  Non0QuadIndex(0, 0, k, k*5, &i0, &i1);
  ASSERT_EQ(0, i0); ASSERT_EQ(2*k, i1);

  Non0QuadIndex(0, 1, k, k*5, &i0, &i1);
  ASSERT_EQ(k, i0); ASSERT_EQ(4, i1);

  Non0QuadIndex(1, 0, k, k*5, &i0, &i1);
  ASSERT_EQ(2, i0); ASSERT_EQ(4, i1);

  Non0QuadIndex(2, 3, k, k*5, &i0, &i1);
  ASSERT_EQ(6, i0); ASSERT_EQ(8, i1);

  Non0QuadIndex(3, 3, k, k*5, &i0, &i1);
  ASSERT_EQ(2*3, i0); ASSERT_EQ(10, i1);

  return 0;
}
int testBSplineSetBasic() {
  PrintTimeStamp(PETSC_COMM_SELF, "set basics", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  
  int order = 3;
  BSS bss; 
  BSSCreate(comm, &bss); 
  BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  ASSERT_EQ(order, bss->order);
  ASSERT_EQ(5, bss->num_ele);
  ASSERT_EQ(5, bss->num_basis);

  ASSERT_EQ(1, bss->b_idx_list[0]);
  ASSERT_EQ(2, bss->b_idx_list[1]);
  ASSERT_EQ(4, bss->b_idx_list[3]);

  PetscReal *zs; BPSGetZs(bps, &zs, NULL);
  for(int i = 0; i < 6; i++)
    ASSERT_DOUBLE_EQ(1.0*i, zs[i]);
  PetscFree(zs);
  ASSERT_DOUBLE_EQ(0.0, bss->ts_r[0]);
  ASSERT_DOUBLE_EQ(0.0, bss->ts_r[1]);
  ASSERT_DOUBLE_EQ(0.0, bss->ts_r[2]);
  ASSERT_DOUBLE_EQ(1.0, bss->ts_r[3]);
  ASSERT_DOUBLE_EQ(2.0, bss->ts_r[4]);
  ASSERT_DOUBLE_EQ(3.0, bss->ts_r[5]);
  ASSERT_DOUBLE_EQ(4.0, bss->ts_r[6]);
  ASSERT_DOUBLE_EQ(5.0, bss->ts_r[7]);
  ASSERT_DOUBLE_EQ(5.0, bss->ts_r[8]);
  ASSERT_DOUBLE_EQ(5.0, bss->ts_r[9]);

  double x = 0.34;
  PetscScalar y1; BSSBasisPsi(bss, 2-1, x, &y1);
  ASSERT_DOUBLE_EQ(0.5*x*x, PetscRealPart(y1));

  ASSERT_DOUBLE_NEAR(0.11270167, bss->xs[0], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(0.5,        bss->xs[1], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(0.88729833, bss->xs[2], pow(10.0, -8.0));

  ASSERT_DOUBLE_NEAR(0.2777777777777, bss->ws[0], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(0.4444444444444, bss->ws[1], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(0.2777777777777, bss->ws[14], pow(10.0, -8.0));

  ASSERT_DOUBLE_NEAR(0.20635083, bss->vals[0], pow(10.0, -8.0));
  ASSERT_DOUBLE_EQ(0.625, bss->vals[1]);
  ASSERT_DOUBLE_EQ(0.0,   bss->vals[10]);
  ASSERT_DOUBLE_EQ(0.125, bss->vals[15*2+4]);

  ASSERT_DOUBLE_NEAR(1.66189500386, bss->derivs[0], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(-0.5         , bss->derivs[4], pow(10.0, -8.0));
  ASSERT_DOUBLE_NEAR(-0.887298334621, bss->derivs[21],pow(10.0, -8.0));

  BSSDestroy(&bss);
  return 0;
}
int testBSplinePsi() {

  PetscErrorCode ierr;
  PrintTimeStamp(PETSC_COMM_SELF, "Psi", NULL);
  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  int n; BSSGetSize(bss, &n);
  if(n <= 0){
    SETERRQ(comm, 1, "n is 0 or negative");
  }
  Vec c;
  ierr = VecCreate(comm, &c); CHKERRQ(ierr);
  ierr = VecSetSizes(c, PETSC_DECIDE, n);CHKERRQ(ierr);
  ierr = VecSetUp(c); CHKERRQ(ierr);

  PetscScalar *ptr_c;
  ierr = VecGetArray(c, &ptr_c); CHKERRQ(ierr);
  for(int i = 0; i < n; i++) {
    ptr_c[i] = i + 0.1;
  }
  ierr = VecRestoreArray(c, &ptr_c); CHKERRQ(ierr);

  PetscReal x = 1.1;
  

  Vec xs; VecCreate(comm, &xs); VecSetSizes(xs, PETSC_DECIDE, 1);
  ierr = VecSetUp(xs); CHKERRQ(ierr);
  ierr = VecSetValue(xs, 0, x, INSERT_VALUES); CHKERRQ(ierr);
  //  Vec xs; VecCreate(comm, &xs); VecSetSizes(xs, PETSC_DEFAULT, 1);
  //  VecSetUp(xs);
  //  ierr = VecSetValue(xs, 0, x, INSERT_VALUES); CHKERRQ(ierr);

  Vec ys;
  ierr = VecDuplicate(xs, &ys); CHKERRQ(ierr);
  ierr = BSSPsi(bss, c, xs, ys); CHKERRQ(ierr);

  Vec dys;
  ierr = VecDuplicate(xs, &dys); CHKERRQ(ierr);
  ierr = BSSDerivPsi(bss, c, xs, dys); CHKERRQ(ierr);

  PetscScalar y;
  ierr = BSSPsiOne(bss, c, x, &y); CHKERRQ(ierr);

  PetscScalar *ptr_y;
  ierr = VecGetArray(ys, &ptr_y); CHKERRQ(ierr);
  ASSERT_SCALAR_EQ(y, ptr_y[0]);
  ierr = VecRestoreArray(ys, &ptr_y); CHKERRQ(ierr);

  PetscScalar dy;
  ierr = BSSDerivPsiOne(bss, c, x, &dy); CHKERRQ(ierr);

  PetscScalar *ptr_dy;
  ierr = VecGetArray(dys, &ptr_dy); CHKERRQ(ierr);
  ASSERT_SCALAR_EQ(dy, ptr_dy[0]);
  ierr = VecRestoreArray(dys, &ptr_dy); CHKERRQ(ierr);

  BSSDestroy(&bss);
  VecRestoreArray(c, &ptr_c);
  VecDestroy(&c);
  VecDestroy(&xs);
  VecRestoreArray(ys, &ptr_y);
  VecDestroy(&ys);
  VecRestoreArray(dys, &ptr_dy);
  VecDestroy(&dys);
  
  return 0;
}
int testBSplineECS() {

  MPI_Comm comm = PETSC_COMM_SELF;
  int order = 8;

  CScaling scaler;
  CScalingCreate(comm, &scaler); CScalingSetSharpECS(scaler, 50.0, 40.0*M_PI/180.0);
  BPS bps; 
  BPSCreate(comm, &bps); BPSSetLine(bps, 70.0, 71);
  BSS bss;
  BSSCreate(comm, &bss);  BSSSetKnots(bss, order, bps); BSSSetCScaling(bss, scaler);

  BSSSetUp(bss);

  PetscScalar y;
  BSSBasisPsi(bss, 60, 58.0, &y);
  ASSERT_DOUBLE_NEAR(0.479365, creal(y), 0.00001);
  ASSERT_DOUBLE_EQ(0.0, cimag(y));

  BSSDestroy(&bss);

  return 0;
}
int testBSplineSetSR1Mat() {
  PrintTimeStamp(PETSC_COMM_SELF, "S", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  // compute S matrix
  Mat S;
  PetscErrorCode ierr;
  ierr = BSSCreateR1Mat(bss, &S);
  ierr = BSSSR1Mat(bss, S);  CHKERRQ(ierr);

  // Get structure and check
  int m, n;
  ierr = MatGetSize(S, &m, &n); CHKERRQ(ierr);
  ASSERT_EQ(5, m);
  ASSERT_EQ(5, n);

  // value check
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  ierr =   MatGetRow(S, 0, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(3, ncols);
  ASSERT_EQ(0, cols[0]);
  ASSERT_DOUBLE_EQ(1.0/3.0, row[0]);
  ASSERT_EQ(1, cols[1]);
  ASSERT_DOUBLE_NEAR(0.2083333333, row[1], pow(10.0, -10.0));
  ASSERT_EQ(2, cols[2]);
  ASSERT_DOUBLE_NEAR(0.0083333333, row[2], pow(10.0, -10.0));
  ierr = MatRestoreRow(S, 0, &ncols, &cols, &row); CHKERRQ(ierr);

  ierr = MatGetRow(S, 2, &ncols, &cols, &row); CHKERRQ(ierr);
  ASSERT_EQ(5, ncols);
  for(int i = 0; i < 5; i++)
    ASSERT_EQ(i, cols[i]);
  ASSERT_DOUBLE_NEAR(0.0083333333333, row[0], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.2166666666666, row[1], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.55, row[2], pow(10.0, -10.0));  
  ASSERT_DOUBLE_NEAR(0.2166666666666, row[3], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.0083333333, row[4], pow(10.0, -10.0));
  ierr = MatRestoreRow(S, 3, &ncols, &cols, &row); CHKERRQ(ierr);

  // Finalize
  MatDestroy(&S);
  BSSDestroy(&bss);

  return 0;
}
int testBSplineSetD2R1Mat() {
  PrintTimeStamp(PETSC_COMM_SELF, "D2", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  // compute matrix
  Mat M; BSSCreateR1Mat(bss, &M); BSSD2R1Mat(bss, M);

  // value check
  const PetscScalar *row;
  MatGetRow(M, 1, NULL, NULL, &row);
  ASSERT_DOUBLE_NEAR(0.1666666666, row[0], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(-1.0, row[1], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.3333333333, row[2], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.16666666666, row[3], pow(10.0, -10.0));
  MatRestoreRow(M, 1, NULL, NULL, &row);

  MatGetRow(M, 2, NULL, NULL, &row);
  ASSERT_DOUBLE_NEAR(0.1666666666666, row[0], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.3333333333333, row[1], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(-1.0, row[2] , pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.3333333333333, row[3], pow(10.0, -10.0));
  ASSERT_DOUBLE_NEAR(0.1666666666666, row[4], pow(10.0, -10.0));
  MatRestoreRow(M, 2, NULL, NULL, &row);

  // Finalize
  MatDestroy(&M);
  BSSDestroy(&bss);

  return 0;
}
int testBSplineSetENMatR1Mat() {
  PrintTimeStamp(PETSC_COMM_SELF, "EN", NULL);
 
  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  // compute matrix
  Mat M; BSSCreateR1Mat(bss, &M); BSSENR1Mat(bss, 2, 0.7, M);

  // value check
  const PetscScalar *row;
  MatGetRow(M, 4, NULL, NULL, &row);
  ASSERT_DOUBLE_NEAR(0.0000964408077, row[0], pow(10.0, -8));
  ASSERT_DOUBLE_NEAR(0.00168615366, row[1], pow(10.0, -8));
  ASSERT_DOUBLE_NEAR(0.00211136230, row[2], pow(10.0, -8));
  MatRestoreRow(M, 4, NULL, NULL, &row);

  // Finalize
  MatDestroy(&M);
  BSSDestroy(&bss);

  return 0;
}
int testBSplineSetEE() {
  PrintTimeStamp(PETSC_COMM_SELF, "EE", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);
  
  Mat ee; BSSCreateR2Mat(bss, &ee); BSSEER2Mat(bss, 0, ee);

  // size check
  PetscInt n, m;
  MatGetSize(ee, &n, &m);
  ASSERT_EQ(25, n); ASSERT_EQ(25, m);

  // value check
  const PetscScalar *row;
  PetscInt ncols;
  const PetscInt *cols;
  MatGetRow(ee, 1, &ncols, &cols, &row);
  ASSERT_DOUBLE_NEAR(0.0709681582, row[0], pow(10.0, -8));
  ASSERT_DOUBLE_NEAR(0.129990244, row[1], pow(10.0, -8));
  ASSERT_DOUBLE_NEAR(0.0371913912, row[2], pow(10.0, -8));
  MatRestoreRow(ee, 1, NULL, NULL, &row);

  /*
    [7.09681582e-02,   1.29990244e-01,   3.71913912e-02,
    1.11566485e-03,   0.00000000e+00,   3.85977593e-02,
    7.84793730e-02,   2.32290473e-02,   6.97290528e-04,
    0.00000000e+00,   1.16319501e-03,   2.84290410e-03,
    9.23893753e-04,   2.78916211e-05,   0.00000000e+00,
    0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    0.00000000e+00]
  */

  // Finalize
  MatDestroy(&ee);
  BSSDestroy(&bss);
  
  return 0;
}
int testBSplineSetEE_time() {
  PrintTimeStamp(PETSC_COMM_SELF, "EE time", NULL);

  time_t t0, t1;
  MPI_Comm comm = PETSC_COMM_SELF;

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 10.0, 31);
  int order = 4;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  t0 = clock();
  Mat ee; BSSCreateR2Mat(bss, &ee); BSSEER2Mat(bss, 0, ee);
  t1 = clock();

  if(getenv("SHOW_DEBUG"))
    PetscPrintf(PETSC_COMM_SELF, "t_new = %f\n", ((double)(t1-t0)/CLOCKS_PER_SEC));

  MatDestroy(&ee);
  BSSDestroy(&bss);
  return 0;
}
int testBSplineSetNE_time() {
  PrintTimeStamp(PETSC_COMM_SELF, "EN time", NULL);

PrintTimeStamp(PETSC_COMM_SELF, "EE time", NULL);

  time_t t0, t1;
  MPI_Comm comm = PETSC_COMM_SELF;

  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 200);
  int order = 10;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetUp(bss);

  t0 = clock();
  Mat en; BSSCreateR1Mat(bss, &en); BSSENR1Mat(bss, 3, 1.3, en);
  t1 = clock();

  if(getenv("SHOW_DEBUG"))
    PetscPrintf(comm, "t(n-e) = %f\n", ((double)(t1-t0)/CLOCKS_PER_SEC));

  MatDestroy(&en);
  BSSDestroy(&bss);
  return 0;
}
int testBSplineHAtom() {
  PrintTimeStamp(PETSC_COMM_SELF, "H atom", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetExp(bps, 30.0, 61, 3.0);
  int order = 5;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);

  Mat H; BSSCreateR1Mat(bss, &H);
  Mat S; BSSCreateR1Mat(bss, &S);
  Mat V; BSSCreateR1Mat(bss, &V);

  BSSD2R1Mat(bss, H);
  MatScale(H, -0.5);

  BSSENR1Mat(bss, 0, 0.0, V);
  MatAXPY(H, -1.0, V, DIFFERENT_NONZERO_PATTERN);

  BSSSR1Mat(bss, S);

  // -- initial space --
  Pot psi0; PotCreate(comm, &psi0); PotSetSlater(psi0, 2.0, 1, 1.1);
  int n_init_space = 1;
  Vec *xs; PetscMalloc1(n_init_space, &xs);
  MatCreateVecs(H, &xs[0], NULL);
  BSSPotR1Vec(bss, psi0, xs[0]);

  EEPS eps; EEPSCreate(comm, &eps);
  EEPSSetOperators(eps, H, S);
  //  EPSSetType(eps->eps, EPSJD);
  EPSSetInitialSpace(eps->eps, 1, xs);
  EEPSSetTarget(eps, -0.6); 

  //  EPSSetInitialSpace(eps->eps, 1, xs);
  
  EEPSSolve(eps);

  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps->eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  EPSGetEigenpair(eps->eps, 0, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.5,  kr, pow(10.0, -6.0));

  Vec cs;
  MatCreateVecs(H, &cs, NULL);
  EEPSGetEigenvector(eps, 0, cs);
  PetscReal x=1.1;
  PetscScalar y=0.0;
  PetscScalar dy=0.0;
  BSSPsiOne(bss, cs, x, &y);
  BSSDerivPsiOne(bss, cs, x, &dy);
  ASSERT_DOUBLE_NEAR(creal(y), 2.0*x*exp(-x), pow(10.0, -6));
  ASSERT_DOUBLE_NEAR(creal(dy), 2.0*exp(-x)-2.0*x*exp(-x), pow(10.0, -6));

  VecDestroy(&xs[0]);
  PetscFree(xs);
  PFDestroy(&psi0);
  BSSDestroy(&bss);
  MatDestroy(&H);
  MatDestroy(&V);
  MatDestroy(&S);
  EEPSDestroy(&eps);
  VecDestroy(&cs);
  
  return 0;
}
int testBSplinePot() {
  PrintTimeStamp(PETSC_COMM_SELF, "pot", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
  Pot pot; PotCreate(comm, &pot); PotSetPower(pot, 1.0, -2.0);

  Mat V; BSSCreateR1Mat(bss, &V);
  Mat U; BSSCreateR1Mat(bss, &U);
  BSSR2invR1Mat(bss, V);
  BSSPotR1Mat(bss, pot, U);
  
  MatAXPY(V, -1.0, U, SAME_NONZERO_PATTERN);
  PetscReal v;
  MatNorm(V, NORM_1, &v);
  ASSERT_DOUBLE_EQ(0.0, v);
  
  BSSDestroy(&bss);
  PFDestroy(&pot);
  MatDestroy(&V);
  MatDestroy(&U);
  return 0;
}
int testBSplinePot2() {
  PrintTimeStamp(PETSC_COMM_SELF, "pot2", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 8);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
  Pot pot; PotCreate(comm, &pot); PotSetCoulombNE(pot, 2, 1.5, 1.0);

  //  POTView(pot);
  Mat V; BSSCreateR1Mat(bss, &V);
  Mat U; BSSCreateR1Mat(bss, &U);
  BSSENR1Mat(bss, 2, 1.5, V);
  BSSPotR1Mat(bss, pot, U);

  MatAXPY(V, -1.0, U, SAME_NONZERO_PATTERN);
  PetscReal v;
  MatNorm(V, NORM_1, &v);
  ASSERT_DOUBLE_EQ(0.0, v);

  BSSDestroy(&bss);
  MatDestroy(&V);
  MatDestroy(&U);
  PFDestroy(&pot);
  
  return 0;
}
int testSlaterPotWithECS() {
  PrintTimeStamp(PETSC_COMM_SELF, "ECS", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 101);
  CScaling scaler; CScalingCreate(comm, &scaler); 
  CScalingSetSharpECS(scaler, 60.0, 20.0*M_PI/180.0);

  int order = 5;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetCScaling(bss, scaler);   BSSSetUp(bss);
  Pot slater; PotCreate(comm, &slater); PotSetSlater(slater, 7.5, 2, 1.0);

  if(getenv("SHOW_DEBUG"))
    BSSView(bss, PETSC_VIEWER_STDOUT_SELF);

  Mat H; BSSCreateR1Mat(bss, &H); 
  Mat V; BSSCreateR1Mat(bss, &V); BSSPotR1Mat(bss, slater, V);
  Mat S; BSSCreateR1Mat(bss, &S); BSSSR1Mat(bss, S);

  BSSD2R1Mat(bss, H);
  MatScale(H, -0.5);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);

  EEPS eps; EEPSCreate(comm, &eps);
  EEPSSetOperators(eps, H, S);
  EEPSSetTarget(eps, 3.4);
  EPSSetDimensions(eps->eps, 10, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSSetTolerances(eps->eps, PETSC_DEFAULT, 1000);
  //  EPSSetType(eps, EPSARNOLDI);

  EEPSSolve(eps);

  PetscInt nconv;
  PetscScalar kr;
  EPSGetConverged(eps->eps, &nconv);
  
  ASSERT_TRUE(nconv > 0);
  if(getenv("SHOW_DEBUG"))
    for(int i = 0; i < nconv; i++) {
      EPSGetEigenpair(eps->eps, i, &kr, NULL, NULL, NULL);
      PetscPrintf(comm, "%f, %f\n", PetscRealPart(kr), PetscImaginaryPart(kr));
    }

  EPSGetEigenpair(eps->eps, 0, &kr, NULL, NULL, NULL);

  PFDestroy(&slater); BSSDestroy(&bss);  EEPSDestroy(&eps);
  MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  
  //  ASSERT_DOUBLE_NEAR(-0.0127745, PetscImaginaryPart(kr), pow(10.0, -4.0));
  //  ASSERT_DOUBLE_NEAR(3.4263903, PetscRealPart(kr), pow(10.0, -4.0));  
  return 0;
}
int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);

  testNumBSpline();
  testCalcBSpline();
  testCalcBSplineECS();
  testNon0QuadIndex();
  testNon0QuadIndex2();
  //  testBSplineECS();

  testBSplineSetBasic();
  testBSplinePsi();
  testBSplineSetSR1Mat();
  testBSplineSetD2R1Mat();
  testBSplineSetENMatR1Mat();

  testBSplineSetEE();

  testBSplineSetEE_time();
  testBSplineSetNE_time();

  testBSplineHAtom();

  testBSplinePot();
  testBSplinePot2();

  testSlaterPotWithECS();
  
  SlepcFinalize();
  return 0;
}


