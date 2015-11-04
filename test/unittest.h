#ifndef UNITTEST_H
#define UNITTEST_H

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

#define ASSERT_NON_NULL(_a)		\
    {\
      void* a = (_a);	\
      if (a == NULL) {				\
	printf("%s:%d: NULL Error: \n", __FILE__, __LINE__);	\
	return 0xffffffffu;	    \
	}\
    }


#endif
