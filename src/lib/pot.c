#include <rescol/pot.h>

typedef struct {
  PetscScalar a;
} Harm;
PetscErrorCode HarmApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  Harm *harm = (Harm*)ctx;
  for(int i = 0; i < n; i++)
    y[i] = harm->a * 0.5 * x[i] * x[i];
  return 0;
}
PetscErrorCode HarmView(void *ctx, PetscViewer v) {
  Harm *self = (Harm*)ctx;
  PetscViewerASCIIPrintf(v, "<<< POT(harm) <<<\n");  
  PetscViewerASCIIPrintf(v, "       a  2 \n");
  PetscViewerASCIIPrintf(v, "v(x) = - x  \n");
  PetscViewerASCIIPrintf(v, "       2    \n");
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "<<< POT(harm) <<<\n");  
  return 0;
}
PetscErrorCode HarmDestroy(void *ctx) {
  PetscErrorCode ierr;
  Harm *self = (Harm*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

typedef struct {
  PetscScalar a;
  PetscInt n;
} Power;
PetscErrorCode PowerApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  Power *self = (Power*)ctx;
  for(int i = 0; i < n; i++)
    y[i] = self->a * pow(x[i], self->n);
  return 0;
}
PetscErrorCode PowerView(void *ctx, PetscViewer v) {
  Power *self = (Power*)ctx;
  PetscViewerASCIIPrintf(v, ">>> POT(Power) >>>\n");
  PetscViewerASCIIPrintf(v, "          n \n");
  PetscViewerASCIIPrintf(v, "v(x) = A x  \n");
  PetscViewerASCIIPrintf(v, "            \n");
  PetscViewerASCIIPrintf(v, "A = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "n = %d\n", self->n);
  PetscViewerASCIIPrintf(v, "<<< POT(Power) <<<\n");  
  return 0;
}
PetscErrorCode PowerDestroy(void *ctx) {
  PetscErrorCode ierr;
  Power *self = (Power*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

typedef struct {
  PetscInt q;
  PetscScalar a;
  PetscScalar zz;
} CoulombNE;
PetscErrorCode CoulombNEApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  CoulombNE *self = (CoulombNE*)ctx;
  PetscScalar a = self->a;
  PetscReal ra = PetscRealPart(a);
  for(int i = 0; i < n; i++) {
    PetscScalar s, g;
    if(ra < PetscRealPart(x[i])) {
      s = a; 
      g = x[i];
    } else {
      s = x[i];
      g = a;
    }
    y[i] = self->zz * pow(s/g, self->q) / g;
  }
  return 0;
}
PetscErrorCode CoulombNEView(void *ctx, PetscViewer v) {
  CoulombNE *self = (CoulombNE*)ctx;
  PetscViewerASCIIPrintf(v, ">>> POT(CoulombNE) >>>\n");
  PetscViewerASCIIPrintf(v, "           s^q   \n");
  PetscViewerASCIIPrintf(v, "v(x) = z ---------\n");
  PetscViewerASCIIPrintf(v, "          g^{q+1} \n");
  PetscViewerASCIIPrintf(v, "s = min{a, x}\n");
  PetscViewerASCIIPrintf(v, "g = max{a, x}\n");
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "q = %d\n", self->q);
  PetscViewerASCIIPrintf(v, "z = %f\n", self->zz);
  PetscViewerASCIIPrintf(v, "<<< POT(CoulombNE) <<<\n");  
  return 0;
}
PetscErrorCode CoulombNEDestroy(void *ctx) {
  PetscErrorCode ierr;
  CoulombNE *self = (CoulombNE*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

typedef struct {
  PetscScalar a;
  PetscInt n;
  PetscScalar z;
} Slater;
PetscErrorCode SlaterApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  Slater *self = (Slater*)ctx;
  for(int i = 0; i < n; i++) {
    y[i] = self->a * pow(x[i], self->n) * exp(-self->z*x[i]);
  }
  return 0;
}
PetscErrorCode SlaterView(void *ctx, PetscViewer v) {
  Slater *self = (Slater*)ctx;
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  PetscViewerType type;     PetscViewerGetType(v, &type);
  PetscViewerFormat format; PetscViewerGetFormat(v, &format);

  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);    

  PetscViewerASCIIPrintf(v, "type: Potential(Slater) \n");
  PetscViewerASCIIPrintf(v, "         n -zx \n");
  PetscViewerASCIIPrintf(v, "v(x) = Ax e    \n");
  PetscViewerASCIIPrintf(v, "A = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "n = %d\n", self->n);
  PetscViewerASCIIPrintf(v, "z = %f\n", self->z);

  return 0;
}
PetscErrorCode SlaterDestroy(void *ctx) {
  PetscErrorCode ierr;
  Slater *self = (Slater*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

typedef struct {
  PetscScalar D0;
  PetscScalar a;
  PetscScalar Re;
} Morse;
PetscErrorCode MorseApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  Morse *self = (Morse*)ctx;
  for(int i = 0; i < n; i++) {
    PetscScalar t = 1.0 - exp(-self->a*(x[i]-self->Re));
    y[i] = self->D0 * t * t;
  }
  return 0;
}
PetscErrorCode MorseView(void *ctx, PetscViewer v) {
  Morse *self = (Morse*)ctx;
  PetscViewerASCIIPrintf(v, ">>> POT(Morse) >>>\n");
  PetscViewerASCIIPrintf(v, "                           2\n");
  PetscViewerASCIIPrintf(v, "v(x) = D0 (1-exp[-a(x-Re)]) \n");
  PetscViewerASCIIPrintf(v, "D0 = %f\n", self->D0);
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "Re = %f\n", self->Re);
  PetscViewerASCIIPrintf(v, "<<< POT(Morse) <<<\n");  
  return 0;
}
PetscErrorCode MorseDestroy(void *ctx) {
  PetscErrorCode ierr;
  Morse *self = (Morse*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PotCreate(MPI_Comm comm, Pot *p_self) {
  PetscErrorCode ierr;
  ierr = PFCreate(comm, 1, 1, p_self); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode PotSetHarm(Pot self, PetscScalar a) {
  Harm *ctx; PetscNew(&ctx);
  ctx->a = a;
  PFSet(self, HarmApply, NULL, HarmView, HarmDestroy, ctx);
  return 0;
}
PetscErrorCode PotSetPower(Pot self, PetscScalar a, PetscInt n) {
  Power *ctx; PetscNew(&ctx);
  ctx->a = a; ctx->n = n;
  PFSet(self, PowerApply, NULL, PowerView, PowerDestroy, ctx);
  return 0;
}
PetscErrorCode PotSetCoulombNE(Pot self, int q, PetscScalar a, PetscScalar zz) {
  CoulombNE *ctx; PetscNew(&ctx);
  ctx->q = q; ctx->a = a; ctx->zz=zz;
  PFSet(self, CoulombNEApply, NULL, CoulombNEView, CoulombNEDestroy, ctx);
  return 0;
}
PetscErrorCode PotSetSlater(Pot self, PetscScalar a, int n, PetscScalar z) {
  Slater *ctx; PetscNew(&ctx);
  ctx->a = a; ctx->n = n; ctx->z = z;
  PFSet(self, SlaterApply, NULL, SlaterView, SlaterDestroy, ctx);
  return 0;
}
PetscErrorCode PotSetMorse(Pot self, PetscScalar D0, PetscScalar a, PetscScalar Re){
  Morse *ctx; PetscNew(&ctx);
  ctx->D0 = D0; ctx->a = a; ctx->Re = Re;
  PFSet(self, MorseApply, NULL, MorseView, MorseDestroy, ctx);
  return 0;
}
PetscErrorCode PotSetFromOptions(Pot self) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  
  char type[10];
  PetscErrorCode ierr;
  PetscBool find;

  ierr = PetscOptionsGetString(NULL, "-pot_type", type, 10, &find); CHKERRQ(ierr);

  if(!find) 
    SETERRQ(comm, 1, "options -pot_type is not found");
  
  if(strcmp(type, "harm") == 0) {
    PetscReal a;
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_a is not found");
    PotSetHarm(self, a);

  } else if(strcmp(type, "slater") == 0) {
    PetscReal v0, z;
    PetscInt n;
    ierr = PetscOptionsGetReal(NULL, "-pot_v0", &v0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_v0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_z", &z, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_z is not found");
    ierr = PetscOptionsGetInt(NULL, "-pot_n", &n, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_n is not found");
    PotSetSlater(self, v0, n, z);

  } else if(strcmp(type, "morse") == 0) {
    PetscReal D0, a, Re;
    ierr = PetscOptionsGetReal(NULL, "-pot_D0", &D0, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_D0 is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_a", &a, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_z is not found");
    ierr = PetscOptionsGetReal(NULL, "-pot_Re", &Re, &find); CHKERRQ(ierr);
    if(!find)
      SETERRQ(comm, 1, "-pot_Re is not found");

    PotSetMorse(self, D0, a, Re);

  } else {
    char msg[100]; sprintf(msg, "unsupported pot_type: %s", type);
    SETERRQ(comm, 1, msg);
  }
  return 0;
}

