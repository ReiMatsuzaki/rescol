#ifndef H_OP_H
#define H_OP_H
#ifdef __cplusplus
extern "C" {
#endif 

#include <petscsys.h>
#include <petscviewer.h>
#include <petscpf.h>

/*
  Operator in QM
*/
#define OpD2 "d2"
#define OpPF "pf"
#define OpPartialVee "partial_vee"
typedef char OpType[256];

struct _p_Op {
  MPI_Comm comm;
  OpType type;
  void *ctx;
  PetscErrorCode (*View)(void*, PetscViewer v);
  PetscErrorCode (*Destroy)(void*);  
  
};
typedef struct _p_Op* Op;
  
PetscErrorCode OpCreate(MPI_Comm comm, Op *p_self);
PetscErrorCode OpDestroy(Op *p_self);
  
PetscErrorCode OpView(Op self, PetscViewer v);
PetscErrorCode OpCheckState(Op self);

PetscErrorCode OpSetD2(Op self);
PetscErrorCode OpSetPF(Op self, PF pf);
PetscErrorCode OpSetPartialVee(Op self, int q);
// PetscErrorCode OpGet(Op self, PetscReal *x);

PetscBool OpIsType(Op self, const OpType type);
PetscErrorCode OpGetPF(Op self, PF *pf);
  
#ifdef __cplusplus
}
#endif
#endif
