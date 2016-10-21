#ifndef RANGE_H
#define RANGE_H
#ifdef __cplusplus
extern "C" {
#endif

#include <petscmat.h>

typedef struct {
  MPI_Comm comm;
  PetscReal x0;
  PetscReal x1;
  PetscReal dx;
  PetscInt num;
  char name[100];
  PetscInt index;  
} p_Range;
typedef p_Range* Range;

PetscErrorCode RangeCreate(MPI_Comm comm, Range *p_self);
PetscErrorCode RangeDestroy(Range *p_self);
PetscErrorCode RangeSetFromOptions(Range self, const char prefix[]);
PetscErrorCode RangeView(Range self, PetscViewer v);
PetscErrorCode RangeSet(Range self, PetscReal x0,
			PetscReal x1, PetscInt num);
PetscErrorCode RangeSetName(Range self, const char name[]);
PetscErrorCode RangeSetFromStr(Range self, const char str[]);
PetscReal RangeGetVal(Range self, int i );

// ---- loop ----
PetscErrorCode RangeInit(Range self);
PetscBool RangeNext(Range self, PetscReal* x);

#ifdef __cplusplus
}
#endif
#endif
