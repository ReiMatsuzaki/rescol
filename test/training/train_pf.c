#include <petscpf.h>

static char help[] = "solve one particle eigen energy problem";


typedef struct {
  PetscScalar a;
  PetscScalar b;
} HarmonicCtx;

PetscErrorCode HarmonicApply(void* ctx, PetscInt n, 
			  const PetscScalar *xs, PetscScalar *ys) {
  HarmonicCtx *harm =  (HarmonicCtx*)ctx;
  PetscScalar a = harm->a;
  PetscScalar b = harm->b;
  ys[0] = a*xs[0]*xs[0] + b;
  ys[1] = xs[0] * b;
  return 0;
}
PetscErrorCode HarmonicApply2(void* ctx, PetscInt n, 
			      const PetscScalar *xs, PetscScalar *ys) {
  HarmonicCtx *harm =  (HarmonicCtx*)ctx;
  PetscScalar a = harm->a;
  for(int i = 0; i < n; i++)
    ys[i] = a*xs[i]*xs[i];
  return 0;
}
PetscErrorCode HarmonicDestroy(void *ctx) {
  PetscErrorCode ierr;
  HarmonicCtx *self = (HarmonicCtx*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode HarmonicView(void* ctx, PetscViewer v) {
  
  PetscViewerType type;
  PetscViewerGetType(v, &type);
  if(strcmp(type, "ascii") != 0) 
    SETERRQ(PETSC_COMM_SELF, 1, "unsupported type");
  
  HarmonicCtx *harm =  (HarmonicCtx*)ctx;
  PetscScalar a = harm->a;
  PetscScalar b = harm->b;
  PetscViewerASCIIPrintf(v, "a*x*x+b\n");
  PetscViewerASCIIPrintf(v, "a=%f\n", a);
  PetscViewerASCIIPrintf(v, "b=%f\n", b);
  return 0;
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;

  ierr = PetscInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);  

  PF pf;
  HarmonicCtx *harm; PetscNew(&harm);
  HarmonicCtx tmp = {1.3,1.2};
  *harm = tmp;
  PFCreate(comm, 1, 2, &pf);
  PFSet(pf, HarmonicApply, NULL, HarmonicView, HarmonicDestroy, harm);
  
  PetscScalar xs[1] = {1.1};
  PetscScalar ys[2] = {0.1, 0.1};
  PFApply(pf, 1, xs, ys);

  printf("%f\n", PetscRealPart(xs[0]));
  printf("%f, %f\n", PetscRealPart(ys[0]), 1.3*1.1*1.1+1.2);
  printf("%f, %f\n", PetscRealPart(ys[1]), 1.1 * 1.2);
  PFView(pf, PETSC_VIEWER_STDOUT_SELF);  
  ierr = PFDestroy(&pf); CHKERRQ(ierr);

  PF pf2;
  HarmonicCtx *harm2; PetscNew(&harm2);
  harm2->a = 1.1; harm2->b = 1.2;
  PFCreate(comm, 1, 1, &pf2); 
  PFSet(pf2, HarmonicApply2, NULL, NULL, HarmonicDestroy, harm2);

  PetscScalar xs2[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  PetscScalar ys2[6];
  PFApply(pf2, 6, xs2, ys2);
  printf("Train2\n");
  for(int i = 0; i < 6; i++) {
    printf("%f\n", PetscRealPart(ys2[i]));
  }
  ierr = PFDestroy(&pf2); CHKERRQ(ierr);  
  
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
