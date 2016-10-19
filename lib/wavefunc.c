#include <rescol/wavefunc.h>

typedef struct {
  int n;
  int L;
  PetscScalar z;
} HEig;
PetscErrorCode HEigApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  HEig *self = (HEig*)ctx;
  for(int i = 0; i < n; i++) {
    y[i] = 2.0 * pow(self->z, 1.5)* x[i] * exp(-self->z*x[i]);
  }
  return 0;
}
PetscErrorCode HEigView(void *ctx, PetscViewer v) {
  HEig *self = (HEig*)ctx;
  PetscViewerASCIIPrintf(v, "type: WaveFunc(Hydrogen Eigen Function) \n");
  PetscViewerASCIIPrintf(v, "n = %d\n", self->n);
  PetscViewerASCIIPrintf(v, "L = %d\n", self->L);  
  PetscViewerASCIIPrintf(v, "z = %f\n", PetscRealPart(self->z));  
  return 0;
}
PetscErrorCode HEigDestroy(void *ctx) {
  PetscErrorCode ierr;
  HEig *self = (HEig*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode WaveFuncCreate(MPI_Comm comm, WaveFunc *p_self) {
  PFCreate(comm, 1, 1, p_self);
  return 0;
}
PetscErrorCode WaveFuncSetHEig(WaveFunc self, int n, int L, PetscScalar z) {
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  if(n != 1 || L != 0)
    SETERRQ(comm, 1, "not implemented yet");

  HEig *ctx; PetscNew(&ctx);
  ctx->n = n;
  ctx->L = L;
  ctx->z = z;
  PFSet(self, HEigApply, 0, HEigView, HEigDestroy, ctx);
  return 0;
}
PetscErrorCode WaveFuncSetFromOptions(WaveFunc self) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  
  char type[10];
  PetscErrorCode ierr;
  PetscBool find;

  ierr = PetscOptionsGetString(NULL, "-wavefunc_type", type, 10, &find); CHKERRQ(ierr);

  if(strcmp(type, "heig") == 0) {
    PetscReal z = 1.0;
    PetscInt n = 1;
    PetscInt L = 0;
    ierr = PetscOptionsGetReal(NULL, "-wavefunc_heig_z", &z, &find); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, "-wavefunc_heig_n", &n, &find); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, "-wavefunc_heig_l", &L, &find); CHKERRQ(ierr);
    ierr = WaveFuncSetHEig(self, n, L, z);
  } else
    SETERRQ(comm, 1, "unsupported wavefunc type");
  return 0;
}

