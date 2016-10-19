#ifndef R1OP_H
#define R1OP_H
#ifdef __cplusplus
// extern "C" {
#endif 
#include <petscmat.h>
#include <petscpf.h>
#include <rescol/bspline.h>
#include <rescol/dvr.h>
#include <rescol/fd.h>

enum OpR1Type {
  R1OpD2,
  R1OpPF
};
struct _p_OpR1 {
  MPI_Comm comm;
  enum OpR1Type type; 
  void *obj;  
  PetscErrorCode (*CalcBSS)();
  PetscErrorCode (*CalcDVR)();
  PetscErrorCode (*CalcFD)();  
  PetscErrorCode (*View)(void *obj, PetscViewer v);
  PetscErrorCode (*Destroy)(void *obj);
};
typedef struct _p_OpR1* OpR1;
typedef struct {
  int q;
} OpCxtEE;


PetscErrorCode OpR1Create(MPI_Comm comm, OpR1* p_self);
PetscErrorCode OpR1SetD2(OpR1 self);
PetscErrorCode OpR1SetDVR(OpR1 self);
PetscErrorCode OpR1SetFD(OpR1 self);

#ifdef __cplusplus
// }
#endif
#endif
