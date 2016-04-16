#ifndef ANGMOMENT_H
#define ANGMOMENT_H
#ifdef __cplusplus
extern "C" {
#endif 

#include <petscmat.h>
#include <gsl/gsl_sf_coupling.h>

#define SIGMA 0
#define PI 1
#define DELTA 2
#define PHI 3
#define GERADE 11
#define UNGERADE 12
#define PLUS 23
#define MINUS 24
#define ROT_SCALAR 101
#define ROT_VECTOR 102
typedef int RotSym;

// ---- utils ----
PetscReal wigner3j(int a, int b, int c, int d, int f, int e);
PetscReal wigner6j(int a, int b, int c, int d, int f, int e);
PetscBool TriangleQ(int j1, int j2, int j3);
int sign(int n);

#ifdef __cplusplus
}
#endif 
#endif
