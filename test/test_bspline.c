#include <slepceps.h>
#include <time.h>
#include "unittest.h"
#include <rescol/bspline.h>

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

  PetscScalar y = 777.0;
  CalcBSpline(order, ts_r, ts_s, 0, 0.0, &y); ASSERT_DOUBLE_EQ(1.0, y);
  CalcBSpline(order, ts_r, ts_s, 0, 1.0, &y); ASSERT_DOUBLE_EQ(0.0, y);
  CalcBSpline(order, ts_r, ts_s, 6, 4.0, &y); ASSERT_DOUBLE_EQ(0.0, y);
  
  PetscScalar x = 0.34;
  CalcBSpline(order, ts_r, ts_s, 2, x, &y); ASSERT_DOUBLE_EQ(0.5*x*x, y);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, &y); ASSERT_DOUBLE_EQ(x, y);

  x = 2.44;
  CalcBSpline(order, ts_r, ts_s, 2, x, &y); ASSERT_DOUBLE_EQ(0.5*x*x-3*x+4.5, y);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, &y);  ASSERT_DOUBLE_EQ(x-3.0, y);

  x = 3.44;
  CalcBSpline(order, ts_r, ts_s, 2, x, &y); ASSERT_DOUBLE_EQ(0.0, y);
  CalcDerivBSpline(order, ts_r, ts_s, 2, x, &y); ASSERT_DOUBLE_EQ(0.0, y);

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
  ASSERT_DOUBLE_NEAR(-0.887298334621, bss->derivs[21],pow(10.0, -8.0))

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
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 20.0, 21);
  int order = 5;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);

  Mat H; BSSCreateR1Mat(bss, &H);
  Mat S; BSSCreateR1Mat(bss, &S);
  Mat V; BSSCreateR1Mat(bss, &V);

  BSSD2R1Mat(bss, H);
  MatScale(H, -0.5);

  BSSENR1Mat(bss, 0, 0.0, V);
  MatAXPY(H, -1.0, V, SAME_NONZERO_PATTERN);

  BSSSR1Mat(bss, S);

  EPS eps; EPSCreate(comm, &eps);
  EPSSetOperators(eps, H, S); EPSSetProblemType(eps, EPS_GHEP);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);
  EPSSetTarget(eps, -0.6); EPSSetFromOptions(eps);
  EPSSolve(eps);

  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  ASSERT_TRUE(nconv > 0);
  EPSGetEigenpair(eps, 0, &kr, NULL, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.5, kr, pow(10.0, -5.0));

  BSSDestroy(&bss);
  MatDestroy(&H);
  MatDestroy(&V);
  MatDestroy(&S);
  EPSDestroy(&eps);
  return 0;
}
int testBSplinePot() {
  PrintTimeStamp(PETSC_COMM_SELF, "pot", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 6);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
  POT pot; POTCreate(comm, &pot); POTSetPower(pot, 1.0, -2.0);

  Mat V; BSSCreateR1Mat(bss, &V);
  Mat U; BSSCreateR1Mat(bss, &U);
  BSSR2invR1Mat(bss, V);
  BSSPotR1Mat(bss, pot, U);
  
  MatAXPY(V, -1.0, U, SAME_NONZERO_PATTERN);
  PetscReal v;
  MatNorm(V, NORM_1, &v);
  ASSERT_DOUBLE_EQ(0.0, v);
  
  return 0;
}
int testBSplinePot2() {
    PrintTimeStamp(PETSC_COMM_SELF, "pot2", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 5.0, 8);
  int order = 3;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);  BSSSetUp(bss);
  POT pot; POTCreate(comm, &pot); POTSetCoulomb(pot, 2.0, 1.5);

  //  POTView(pot);

  Mat V; BSSCreateR1Mat(bss, &V);
  Mat U; BSSCreateR1Mat(bss, &U);
  BSSENR1Mat(bss, 2, 1.5, V);
  BSSPotR1Mat(bss, pot, U);

  MatAXPY(V, -1.0, U, SAME_NONZERO_PATTERN);
  PetscReal v;
  MatNorm(V, NORM_1, &v);
  ASSERT_DOUBLE_EQ(0.0, v);
  
  return 0;
}
int testSlaterPotWithECS() {
  PrintTimeStamp(PETSC_COMM_SELF, "ECS", NULL);

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 101);
  Scaler scaler; ScalerCreate(comm, &scaler); 
  ScalerSetSharpECS(scaler, 60.0, 20.0*M_PI/180.0);

  int order = 5;
  BSS bss; BSSCreate(comm, &bss); BSSSetKnots(bss, order, bps);
  BSSSetScaler(bss, scaler);   BSSSetUp(bss);
  POT slater; POTCreate(comm, &slater); POTSetSlater(slater, 7.5, 1.0);

  if(getenv("SHOW_DEBUG"))
    BSSView(bss, PETSC_VIEWER_STDOUT_SELF);

  Mat H; BSSCreateR1Mat(bss, &H); 
  Mat V; BSSCreateR1Mat(bss, &V); BSSPotR1Mat(bss, slater, V);
  Mat S; BSSCreateR1Mat(bss, &S); BSSSR1Mat(bss, S);

  BSSD2R1Mat(bss, H);
  MatScale(H, -0.5);
  MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN);

  EPS eps; EPSCreateForBoundState(&eps, comm, H, S, 3.4);
  EPSSetDimensions(eps, 10, PETSC_DEFAULT, PETSC_DEFAULT);
  EPSSetTolerances(eps, PETSC_DEFAULT, 1000);
  //  EPSSetType(eps, EPSARNOLDI);

  EPSSolve(eps);

  PetscInt nconv;
  PetscScalar kr, ki;
  EPSGetConverged(eps, &nconv);
  
  ASSERT_TRUE(nconv > 0);
  if(getenv("SHOW_DEBUG"))
    for(int i = 0; i < nconv; i++) {
      EPSGetEigenpair(eps, i, &kr, &ki, NULL, NULL);
      PetscPrintf(comm, "%f, %f\n", PetscRealPart(kr), PetscImaginaryPart(kr));
    }

  EPSGetEigenpair(eps, 0, &kr, &ki, NULL, NULL);
  ASSERT_DOUBLE_NEAR(-0.0127745, PetscImaginaryPart(kr), pow(10.0, -4.0));
  ASSERT_DOUBLE_NEAR(3.4263903, PetscRealPart(kr), pow(10.0, -4.0));  

  POTDestroy(&slater);
  BSSDestroy(&bss); MatDestroy(&H); MatDestroy(&V); MatDestroy(&S);
  EPSDestroy(&eps);
  return 0;
}
int main(int argc, char **args) {

  SlepcInitialize(&argc, &args, (char*)0, help);

  testNumBSpline();
  testCalcBSpline();
  testNon0QuadIndex();
  testNon0QuadIndex2();
  
  testBSplineSetBasic();
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


