#include <rescol/op.h>

typedef struct {
  PetscReal a;
  int q;
  PetscReal zz;
} OpCtxCoulomb;
PetscErrorCode OpCtxCoulombView(void *ctx, PetscViewer v) {
  PetscErrorCode ierr;
  OpCtxCoulomb *pt = (OpCtxCoulomb*)ctx;
  ierr = PetscViewerASCIIPrintf(v, "OpCoulomb object\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v, "a = %f\n", pt->a); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v, "q = %d\n", pt->q); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(v, "zz= %f\n", pt->zz); CHKERRQ(ierr);  
  return 0;
}
PetscErrorCode OpCtxCoulombDestroy(void *ctx) {
  PetscErrorCode ierr;
  OpCtxCoulomb *pt = (OpCtxCoulomb*)ctx;
  ierr = PetscFree(pt); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode OpCtxPFView(void *ctx, PetscViewer v) {
  PF pf = (PF)ctx;
  PFView(pf, v);
  return 0;
}
PetscErrorCode OpCtxPFDestroy(void *ctx) {
  PF pf = (PF)ctx;
  PFDestroy(&pf);
  return 0;
}

PetscErrorCode OpCreate(MPI_Comm comm, Op *p_self) {
  PetscErrorCode ierr;
  Op self;  
  ierr = PetscNew(&self); CHKERRQ(ierr);
  self->comm = comm;
  *p_self = self;
  return 0;	       
}
PetscErrorCode OpDestroy(Op *p_self) {
  PetscErrorCode ierr;
  Op self = *p_self;

  if(self->Destroy)
    self->Destroy(self->ctx);

  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode OpView(Op self, PetscViewer v) {
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  
  ierr = OpCheckState(self); CHKERRQ(ierr);
  
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);  
  
  if(iascii) {
    PetscViewerASCIIPrintf(v, "Op object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "type: %s\n", self->type);
    if(self->View) 
      ierr = self->View(self->ctx, v);
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {
    
  } else if(isdraw) {
    
  }
  return 0;
  
}
PetscErrorCode OpCheckState(Op self) {
  //  PetscErrorCode ierr;
  if(self == NULL)
    SETERRQ(PETSC_COMM_SELF, 1, "object is null");
  return 0;
}

PetscErrorCode OpSetD2(Op self) {
  //  PetscErrorCode ierr;
  strcpy(self->type, OpD2);
  self->ctx = NULL;
  self->View = NULL;
  self->Destroy = NULL;
  return 0;
}
PetscErrorCode OpSetPF(Op self, PF pf) {
  strcpy(self->type, OpPF);
  self->ctx = pf;
  self->View = OpCtxPFView;
  self->Destroy = OpCtxPFDestroy;
  return 0;
}
PetscErrorCode OpSetPartialVee(Op self, int q) {
  PetscErrorCode ierr;
  OpCtxCoulomb *ctx; 
  ierr = PetscNew(&ctx); CHKERRQ(ierr);

  ctx->q = q;
  ctx->a = 0.0;
  ctx->zz = 1.0;

  strcpy(self->type, OpPartialVee);
  self->ctx = ctx;
  self->View = OpCtxCoulombView;
  self->Destroy = OpCtxCoulombDestroy;
  return 0;
}
//PetscErrorCode OpGet(Op self, PetscReal *x);

PetscBool OpIsType(Op self, const OpType type) {
  return (strcmp(self->type, type) == 0);
}
PetscErrorCode OpGetPF(Op self, PF *pf) {
  if(!OpIsType(self, OpPF))
    SETERRQ(self->comm, 1, "type is not OpPF");

  *pf = (PF)self->ctx;
  return 0;
}



