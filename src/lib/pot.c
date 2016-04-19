#include <string.h>
#include <stdlib.h>
#include <rescol/pot.h>

int str_nele(const char *str, char key) {

  int num = strlen(str);
  int cumsum = 0;
  for(int i = 0; i < num; i++) {
    if(str[i] == key)
      cumsum++;
  }
  return cumsum;
}

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
  PetscViewerASCIIPrintf(v, "type: Potential(Harm)\n");
  PetscViewerASCIIPrintf(v, "       a  2 \n");
  PetscViewerASCIIPrintf(v, "v(x) = - x  \n");
  PetscViewerASCIIPrintf(v, "       2    \n");
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
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
  PetscViewerASCIIPrintf(v, "type: Potential(Power)\n");
  PetscViewerASCIIPrintf(v, "          n \n");
  PetscViewerASCIIPrintf(v, "v(x) = A x  \n");
  PetscViewerASCIIPrintf(v, "            \n");
  PetscViewerASCIIPrintf(v, "A = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "n = %d\n", self->n);
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
  PetscViewerASCIIPrintf(v, "type: Potential(CoulombNE)\n");
  PetscViewerASCIIPrintf(v, "           s^q   \n");
  PetscViewerASCIIPrintf(v, "v(x) = z ---------\n");
  PetscViewerASCIIPrintf(v, "          g^{q+1} \n");
  PetscViewerASCIIPrintf(v, "s = min{a, x}\n");
  PetscViewerASCIIPrintf(v, "g = max{a, x}\n");
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "q = %d\n", self->q);
  PetscViewerASCIIPrintf(v, "z = %f\n", self->zz);
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
  PetscViewerASCIIPrintf(v, "type: Potential(Mourse)\n");
  PetscViewerASCIIPrintf(v, "v(x) = D0 (1-exp[-a(x-Re)]) \n");
  PetscViewerASCIIPrintf(v, "D0 = %f\n", self->D0);
  PetscViewerASCIIPrintf(v, "a = %f\n", self->a);
  PetscViewerASCIIPrintf(v, "Re = %f\n", self->Re);
  return 0;
}
PetscErrorCode MorseDestroy(void *ctx) {
  PetscErrorCode ierr;
  Morse *self = (Morse*)ctx;
  ierr = PetscFree(self); CHKERRQ(ierr);
  return 0;
}

typedef struct {
  int num;
  PF *pfs;
} Combination;
PetscErrorCode CombinationApply(void *ctx, int n, const PetscScalar *x, PetscScalar *y) {
  PetscErrorCode ierr;
  Combination *self = (Combination*)ctx;
  PetscScalar *y_tmp; PetscMalloc1(n, &y_tmp);
  for(int i = 0; i < n; i++)
    y[i] = 0.0;
  for(int j = 0; j < self->num; j++) {
    PFApply(self->pfs[j], n, x, y_tmp);
    for(int i = 0; i < n; i++)
      y[i] += y_tmp[i];
  }
  ierr = PetscFree(y_tmp); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode CombinationView(void *ctx, PetscViewer v) {
  PetscErrorCode ierr;
  Combination *self = (Combination*)ctx;
  PetscViewerASCIIPrintf(v, "type: Combination of PF object\n");
  PetscViewerASCIIPrintf(v, "# of PF: %d\n", self->num);
  for(int i = 0; i < self->num; i++) {
    PetscViewerASCIIPrintf(v, "No %d:\n", i);
    ierr = PFView(self->pfs[i], v); CHKERRQ(ierr);
  }
  return 0;
}
PetscErrorCode CombinationDestroy(void *ctx) {
  PetscErrorCode ierr;
  Combination *self = (Combination*)ctx;
  ierr = PetscFree(self->pfs); CHKERRQ(ierr);
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
PetscErrorCode PotSetCombination(Pot self, int num, Pot *pfs) {
  Combination *ctx; PetscNew(&ctx);
  ctx->num = num;
  PetscMalloc1(num, &ctx->pfs);
  for(int i = 0; i < num; i++)
    ctx->pfs[i] = pfs[i];
  PFSet(self, CombinationApply, NULL, CombinationView, CombinationDestroy, ctx);
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

  } else if(strcmp(type, "coulomb") == 0) {
    PetscReal z = -1.0;
    ierr = PetscOptionsGetReal(NULL, "-pot_z", &z, &find); CHKERRQ(ierr);
    ierr = PotSetCoulombNE(self, 0, 0.0, z); CHKERRQ(ierr);
  } else {
    char msg[100]; sprintf(msg, "unsupported pot_type: %s", type);
    SETERRQ(comm, 1, msg);
  }
  return 0;
}
PetscErrorCode PotSetFromStr_Combination(Pot self, char str[]) {
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  PetscErrorCode ierr;
  int num_prim = str_nele(str, '|') + 1;
  PF *pfs; PetscMalloc1(num_prim, &pfs);

  char **prim;
  PetscMalloc1(num_prim, &prim);
  prim[0] = strtok(str, "|");
  for(int i = 1; i < num_prim; i++) {
    prim[i] = strtok(NULL, "|");
  }

  for(int i = 0; i < num_prim; i++) {

    if(prim[i] == NULL) {
      char msg[100]; 
      sprintf(msg, "null found. i = %d. str = %s", i, str);
      SETERRQ(comm, 1, msg);
    }

    ierr = PotCreate(comm, &pfs[i]); CHKERRQ(ierr);
    ierr = PotSetFromStr(pfs[i], prim[i]); CHKERRQ(ierr);

  }

  ierr = PotSetCombination(self, num_prim, pfs); CHKERRQ(ierr);
  PetscFree(prim);
  return 0;
}
PetscErrorCode PotSetFromStr_Slater(Pot self, char str[]) {
  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  char *err;
  char *str_c;
  double c;
  str_c = strtok(NULL, " ");
  if(str_c == NULL) {
    char msg[100]; sprintf(msg, "c is not found. original string: %s", str);
    SETERRQ(comm, 1, msg);
  }
  c = atof(str_c);
  if(fabs(c) < 0.00000000001)
    SETERRQ(comm, 1, "coefficient c is 0");

  char *str_n;
  int n;
  str_n = strtok(NULL, " ");
  if(str_n == NULL) {
    char msg[100]; sprintf(msg, "n is not found. Original string: %s", str);
    SETERRQ(comm, 1, msg);
  }
  n = strtol(str_n, &err, 10);
  if(err == '\0')
    SETERRQ(comm, 1, "conversion failed for n");
    
  char *str_z;
  double z;
  str_z = strtok(NULL, " ");
  if(str_z == NULL) {
    char msg[100]; sprintf(msg, "z is not found. Original string: %s", str);
    SETERRQ(comm, 1, msg);
  }
  z = atof(str_z);
  if(fabs(z) < 0.00000000001)
    SETERRQ(comm, 1, "orbital exponent z is 0");

  PotSetSlater(self, c, n, z);
  return 0;
}
PetscErrorCode PotSetFromStr_Power(Pot self, char str[]) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  char *err;

  char *str_c;
  double c;
  str_c = strtok(NULL, " ");
  if(str_c == NULL)
    SETERRQ(comm, 1, "c is not found");
  c = atof(str_c);
  if(fabs(c) < 0.00000000001)
    SETERRQ(comm, 1, "coefficient c is 0");
  
  char *str_n;
  int n;
  str_n = strtok(NULL, " ");
  if(str_n == NULL)
    SETERRQ(comm, 1, "n is not found");
  n = strtol(str_n, &err, 10);
  if(err == '\0')
    SETERRQ(comm, 1, "conversion failed for n");
  
  PotSetPower(self, c, n);
  return 0;

}
PetscErrorCode PotSetFromStr(Pot self, const char str_in[]) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  PetscErrorCode ierr;

  char str[101];
  if(strlen(str_in) > 100)
    SETERRQ(comm, 1, "length of str_in must be shorter than 100");
  strcpy(str, str_in);

  // -- check summation --
  if(str_nele(str_in, '|') != 0) {
    ierr = PotSetFromStr_Combination(self, str); CHKERRQ(ierr);
    return 0;
  }

  // -- check mono --
  char *name;
  name = strtok(str, " ");
  if(name == NULL) {
    SETERRQ(comm, 1, "name is not found in str\n"); 
  }

  if(strcmp(name, "sto") == 0) {
    ierr = PotSetFromStr_Slater(self, str); CHKERRQ(ierr);
    return 0;
  }
  else if(strcmp(name, "pow") == 0) {
    ierr = PotSetFromStr_Power(self, str); CHKERRQ(ierr);
    return 0;
  }
  else {
    SETERRQ(comm, 1, "Unsupported name");
  }
}
PetscErrorCode PotSetFromOptions2(Pot self, const char prefix[]) {

  MPI_Comm comm; PetscObjectGetComm((PetscObject)self, &comm);
  PetscErrorCode ierr;
  PetscBool find;

  char option_name[100] = "-";
  strcat(option_name, prefix);
  strcat(option_name, "-pot");

  char option_res[100];

  ierr = PetscOptionsGetString(NULL, option_name, option_res, 100, &find); CHKERRQ(ierr);

  if(!find) {
    char msg[1000];
    sprintf(msg, "Not found option. option_name = %s.", option_name);
    SETERRQ(comm, 1, msg);
  }

  ierr = PotSetFromStr(self, option_res);

  return 0;
}

