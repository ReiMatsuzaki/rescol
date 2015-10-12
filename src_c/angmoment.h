#ifndef ANGMOMENT_H
#define ANGMOMENT_H
#include <petscmat.h>
#include <Python/Python.h>

#define GERADE 1
#define UNGERADE 2

#define GERADE_PLUS 3
#define GERADE_MINUS 4
#define GERADE_ALL 5
#define UNGERADE_PLUS 6
#define UNGERADE_MINUS 7
#define UNGERADE_ALL 8

double w3j(int j1, int j2, int j3, int m1, int m2, int m3);
double w6j(int j1, int j2, int j3, int j4, int j5, int j6);
double PyGaunt(int j1, int j2, int j3, int m1, int m2, int m3);
double Y1RedMat_Yq(int j1, int j2, int j3);
double Y1Mat_Yqk(int j1, int j2, int j3, int m1, int m2, int m3);

struct _p_Y1s {
  int num;
  int* ls;
  int m;
};

typedef struct _p_Y1s* Y1s;

typedef int Y1Sym;
typedef int Y2Sym;

PetscErrorCode Y1sCreate(Y1s *y1s, int l0, int maxl, Y1Sym sym, int m);
PetscErrorCode Y1sDestroy(Y1s *y1s);
PetscErrorCode Y1sView(Y1s y1s);
PetscErrorCode Y1sCreateY1Mat(Y1s ys, MPI_Comm comm, Mat *M);
PetscErrorCode Y1sCalcLambdaY1Mat(Y1s ys, Mat M, InsertMode mode);
PetscErrorCode Y1sCalcPqY1Mat(Y1s ys, int q, Mat M, InsertMode mode);

#endif
