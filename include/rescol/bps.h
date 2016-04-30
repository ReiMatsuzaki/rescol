#ifndef BPS_H
#define BPS_H
#ifdef __cplusplus
extern "C" {
#endif

#include <petscmat.h>

/*
  Break points for finite elements methods
*/

struct _p_BPS {
  char type[10];
  PetscInt num_zs;
  PetscReal *zs;
  MPI_Comm comm;
};

typedef struct _p_BPS* BPS;

// ------- Basic ----------
PetscErrorCode BPSCreate(MPI_Comm comm, BPS *p_self);
PetscErrorCode BPSDestroy(BPS *p_self);
PetscErrorCode BPSView(BPS self, PetscViewer v);
PetscErrorCode BPSCheckState(BPS self);

// ---- Setter ----
PetscErrorCode BPSSetFromOptions(BPS self);
PetscErrorCode BPSSetExp(BPS self, PetscReal zmax, PetscInt num_zs, PetscReal gamma);
PetscErrorCode BPSSetLine(BPS self, PetscReal zmax, PetscInt num_zs);

// ---- Getter ----
PetscErrorCode BPSGetZs(BPS self, PetscReal **zs, PetscInt *num_zs);
PetscErrorCode BPSGetNumEle(BPS self, PetscInt *num_ele);
PetscErrorCode BPSGetZMax(BPS self, PetscReal *zmax);
PetscErrorCode BPSGetEdge(BPS self, int iele, PetscReal *z0, PetscReal *z1);
PetscErrorCode BPSInElementQ(BPS self, int iele, PetscReal x, PetscBool *in_q);

#ifdef __cplusplus
}
#endif
#endif
