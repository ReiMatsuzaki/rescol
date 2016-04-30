#ifndef UNITTEST_H
#define UNITTEST_H

#include <slepceps.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define ASSERT_TRUE(_a) \
    {\
        int assertErrorA = (_a);\
        if (! assertErrorA) {\
	  printf("%s:%d: TRUE Error\n", __FILE__, __LINE__);\
	  printf("expectted: True\n");\
	  printf("actual:    False\n");\
	  return 0xffffffffu;	       \
        }\
    }
#define ASSERT_FALSE(_a) \
    {\
        int assertErrorA = (_a);\
        if (assertErrorA) {\
	  printf("%s:%d: FALSE Error\n", __FILE__, __LINE__);\
	  printf("expectted: False\n");\
	  printf("actual:    True\n");\
	  return 0xffffffffu;	       \
        }\
    }

#define ASSERT_EQ(_a,_b) \
    {\
        int assertErrorA = (_a);\
        int assertErrorB = (_b);\
        if ((assertErrorA) != (assertErrorB)) {\
	  printf("%s:%d: EQ Error: \n", __FILE__, __LINE__);\
	  printf("expectd: %d\n", _a);\
	  printf("actual : %d\n", _b);\
	  return 0xffffffffu;	    \
        }\
    }

#define ASSERT_NOT_EQ(_a,_b) \
    {\
        int assertErrorA = (_a);\
        int assertErrorB = (_b);\
        if ((assertErrorA) == (assertErrorB)) {\
	  printf("%s:%d: NOT_EQ Error: \n", __FILE__, __LINE__); \
	  printf("expectd: %d\n", _a);			     \
	  printf("actual : %d\n", _b);\
	  return 0xffffffffu;	    \
	}\
    }

#define ASSERT_DOUBLE_EQ(_a,_b) \
    {\
      double a = (_a);			\
      double b = (_b);			\
      if (fabs(a-b) > 10.0*DBL_EPSILON) {				\
	printf("%s:%d: DOUBLE_EQ Error: \n", __FILE__, __LINE__);	\
	printf("expectd: %f\n", a);\
	printf("actual : %f\n", b); \
	printf("diff : %e\n", a-b); \
	return 0xffffffffu;	    \
	}\
    }

#define ASSERT_DOUBLE_NEAR(_a,_b, _eps)		\
    {\
      double a = (_a);			\
      double b = (_b);			\
      double eps = (_eps); \
      if (fabs(a-b) > eps) {				\
	printf("%s:%d: DOUBLE_EQ Error: \n", __FILE__, __LINE__);	\
	printf("expectd: %f\n", a);\
	printf("actual : %f\n", b); \
	printf("a-b : %e\n", a-b);  \
	printf("eps : %e\n", eps);  \
	return 0xffffffffu;	    \
	}\
    }

#define ASSERT_SCALAR_EQ(_a,_b) \
    {\
      PetscScalar a = (_a);			\
      PetscScalar b = (_b);			\
      if(cabs(a-b) > 10.0*DBL_EPSILON) {				\
	printf("%s:%d: SCALAR_EQ Error: \n", __FILE__, __LINE__); \
	printf("expectd: %f, %f\n", creal(a), cimag(a));	  \
	printf("actual : %f, %f\n", creal(b), cimag(b));          \
	printf("diff : %e, %e\n", creal(a-b), cimag(a-b));	  \
	return 0xffffffffu;	    \
	}\
    }
#define ASSERT_SCALAR_NEAR(_a,_b,_eps)		\
    {\
      PetscScalar a = (_a);			\
      PetscScalar b = (_b);			\
      PetscReal eps = (_eps);				\
      if(cabs(a-b) > eps) {				\
	printf("%s:%d: SCALAR_NEAR Error: \n", __FILE__, __LINE__); \
	printf("expectd: %f, %f\n", creal(a), cimag(a));	  \
	printf("actual : %f, %f\n", creal(b), cimag(b));          \
	printf("diff : %e, %e\n", creal(a-b), cimag(a-b));	  \
	printf("eps : %e\n", eps);	  \
	return 0xffffffffu;	    \
	}\
    }
PetscErrorCode test_mat_eq(Mat A, Mat B, double eps, const char file[], int line) {
  PetscErrorCode ierr;
  int na, ma, nb, mb;
  MatGetSize(A, &na, &ma);
  MatGetSize(B, &nb, &mb);
  if(na != nb || ma != mb) {
    printf("%s:%d: size mismatch.\n", file, line);
    printf("size of A = (%d, %d)\n", na, ma);
    printf("size of B = (%d, %d)\n", nb, mb);
    return 1;
  }
  if(na != ma ) {
    printf("%s:%d: this test only support square matrix.\n", file, line);
    printf("size of A = (%d, %d)\n", na, ma);
    return 1;
  }
  PetscInt *idx;
  ierr = PetscMalloc1(na, &idx); CHKERRQ(ierr);
  for(int i = 0; i < na; i++) { idx[i] = i; }
  PetscScalar *xa, *xb;
  ierr = PetscMalloc1(na*nb, &xa); CHKERRQ(ierr);
  ierr = MatGetValues(A, na, idx, na, idx, xa); CHKERRQ(ierr);
  
  ierr = PetscMalloc1(na*nb, &xb); CHKERRQ(ierr);
  ierr = MatGetValues(B, na, idx, na, idx, xb); CHKERRQ(ierr);
  for(int i = 0; i < na; i++)
    for(int j = 0; j < na; j++) {      
      PetscScalar va = xa[i+na*j];
      PetscScalar vb = xb[i+na*j];
      if((cabs(va) > 10.0) ?
	 cabs((va-vb)/va) > eps:
	 cabs(va-vb) > eps) {
	printf("%s:%d: element is not equal.\n", file, line);
	printf("(i, j) = (%d, %d)\n", i, j);
	PetscScalar vala = xa[i+na*j];
	PetscScalar valb = xb[i+na*j];
	printf("expectd: %f, %f\n", creal(vala), cimag(vala));
	printf("actual : %f, %f\n", creal(valb), cimag(valb));
	printf("|a-b| : %e\n", cabs(vala-valb));
	printf("eps : %e\n", eps);
	return 1;
      }
    }
  PetscFree(idx);
  PetscFree(xa);
  PetscFree(xb);
  return 0;
}
PetscErrorCode test_vec_eq(Vec A, Vec B, const char file[], int line) {
  int na, nb; VecGetSize(A, &na); VecGetSize(B, &nb);
  if(na != nb) {
    printf("%s:%d: size mismatch.\n", file, line);
    printf("size of A = %d\n", na);
    printf("size of B = %d\n", nb);
    return 1;
  }

  PetscScalar *xa, *xb;
  VecGetArray(A, &xa);
  VecGetArray(B, &xb);
  for(int i = 0; i < na; i++) {
    if(cabs(xa[i]-xb[i]) > 10.0*DBL_EPSILON) {
      printf("%s:%d: element is not near.\n", file, line);
      printf("i = %d\n", i);
      printf("A[i] = %f, %f\n", creal(xa[i]), cimag(xa[i]));
      printf("B[i] = %f, %f\n", creal(xb[i]), cimag(xb[i]));
      return 1;
    }
  }
  VecRestoreArray(A, &xa);
  VecRestoreArray(B, &xb);
  return 0;
}

#define ASSERT_MAT_EQ(_A, _B) \
  {\
    Mat A = (_A); \
    Mat B = (_B); \
    if((!test_mat_eq(A, B, __FILE__, __LINE__))) {		\
      printf("%s:%d: NULL Error: \n", __FILE__, __LINE__);	\
      return 0xffffffffu;					\
    }\
    printf("E\n");\
  }

#define ASSERT_NON_NULL(_a)		\
    {\
      void* a = (_a);	\
      if (a == NULL) {				\
	printf("%s:%d: NULL Error: \n", __FILE__, __LINE__);	\
	return 0xffffffffu;	    \
	}\
    }


#endif
