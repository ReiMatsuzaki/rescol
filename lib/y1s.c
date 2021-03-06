#include <float.h>
#include "../include/math.h"
#include "../include/y1s.h"

// ---- Y1 Functions ----
PetscReal Y1RedYq(int j1, int q, int j2) {
  return sign(j1) * 
    sqrt((2*j1+1)*(2*j2+1)*(2*q+1)/(4*M_PI)) *     
    wigner3j(j1, q, j2, 0, 0, 0);
}
PetscReal Y1EleYqk(int j1, int q, int j2, int m1, int k, int m2) {
  return sign(j1-m1)*
    Y1RedYq(j1, q, j2) * 
    wigner3j(j1, q, j2, -m1, k, m2);
}
PetscReal Y1ElePq(int j1, int q, int j2, int m1, int m2) {
  return Y1EleYqk(j1, q, j2, m1, 0, m2) * sqrt(4.0*M_PI/(2*q+1));
}

// ---- Y1s Method ----
PetscErrorCode Y1sCreate(MPI_Comm comm, Y1s *p_self) {

  PetscErrorCode ierr;
  Y1s self;
  ierr = PetscNew(&self); CHKERRQ(ierr);

  self->comm = comm;
  //  self->prefix = "";
  self->num = 0;
  self->ls = NULL;
  self->m = 0;
  
  *p_self = self;
  return 0;
}
PetscErrorCode Y1sDestroy(Y1s* p_self) {

  PetscFree((*p_self)->ls);
  PetscFree(*p_self);
  return 0;

}

PetscErrorCode Y1sView(Y1s self, PetscViewer v) {

  PetscErrorCode ierr;
  if(self == NULL) {
    SETERRQ(MPI_COMM_SELF, 1, "self is null");
  }
  
  PetscBool iascii, isbinary, isdraw;
  PetscViewerType type;     PetscViewerGetType(v, &type);
  PetscViewerFormat format; PetscViewerGetFormat(v, &format);

  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);    

  if(iascii) {
    PetscViewerASCIIPrintf(v, "Y1s object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "num: %d\n", self->num);
    PetscViewerASCIIPrintf(v, "m: %d\n", self->m);    
    PetscViewerASCIIPrintf(v, "ls: ");
    PetscViewerASCIIUseTabs(v, PETSC_FALSE);
    for(int i = 0; i < self->num; i++) {
      PetscViewerASCIIPrintf(v, "%d ", self->ls[i]);
    }
    PetscViewerASCIIPrintf(v, "\n");
    PetscViewerASCIIUseTabs(v, PETSC_TRUE);
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {
    SETERRQ(self->comm, 1, "not supported");
  } else if(isdraw) {
    SETERRQ(self->comm, 1, "not supported");
  }
  return 0;
}

PetscErrorCode Y1sSet(Y1s self, int m, int g_or_u, int lmax) {

  PetscErrorCode ierr;

  if(g_or_u != GERADE && g_or_u != UNGERADE)
    SETERRQ(self->comm, 1, "illegal g_or_u");

  self->m = m;

  int l0 = (g_or_u == GERADE ? 0 : 1);

  self->num = (lmax-l0)/2 + 1;
  ierr = PetscMalloc1(self->num, &self->ls); CHKERRQ(ierr);  
  
  int i = 0;
  for(int L = l0; L <= lmax; L+=2) {
    self->ls[i] = L; i++;
  }
  return 0;
}
PetscErrorCode Y1sSetOne(Y1s self, int L, int M) {

  PetscErrorCode ierr;
  if(L < 0) 
    SETERRQ(self->comm, 1, "L must be zero or positive integer");

  if(-M < -L || L < M)
    SETERRQ(self->comm, 1, "M must smaller or equal to L");

  self->num = 1;
  ierr = PetscMalloc1(1, &self->ls); CHKERRQ(ierr);
  self->ls[0] = L;
  self->m = M;

  return 0;
}
PetscErrorCode Y1sSetOptionsPrefix(Y1s self, const char prefix[]) {
  strcpy(self->prefix, prefix);
  return 0;
}
PetscErrorCode Y1sSetFromOptions(Y1s self) {

  char rot[10] = "sigma";
  int lmax = 2;
  PetscErrorCode ierr;
  PetscBool find;

  char opt_rot[100];
  sprintf(opt_rot, "-%sy1s_rot", self->prefix);
  PetscOptionsGetString(NULL, NULL, opt_rot, rot, 10, NULL);

  int m = SIGMA;
  if(strcmp(rot, "sigma") == 0)
    m = SIGMA;
  else if(strcmp(rot, "pi") == 0)
    m = PI;
  else if(strcmp(rot, "delta") == 0)
    m = DELTA;
  else if(strcmp(rot, "phi") == 0)
    m = PHI;
  else
    SETERRQ(self->comm, 1, "options -rot <- {sigma, pi, delta, phi}");

  char opt_lmax[100];
  sprintf(opt_lmax, "-%sy1s_lmax", self->prefix);
  PetscOptionsGetInt(NULL, NULL, opt_lmax, &lmax, &find); 

  if(find) {
    
    char parity[10] = "gerade";
    char opt_parity[100];
    int g_or_u = GERADE;

    sprintf(opt_parity, "-%sy1s_parity", self->prefix);
    PetscOptionsGetString(NULL, NULL, opt_parity, parity, 10, &find);
    
    if(strcmp(parity, "gerade") == 0)
      g_or_u = GERADE;
    else if(strcmp(parity, "ungerade") == 0)
      g_or_u = UNGERADE;
    else
      SETERRQ(self->comm, 1, "options -parity <- {gerade, ungerade}");
    if(lmax < 0)
      SETERRQ(self->comm, 1, "options lmax must non negative integer");
    
    ierr = Y1sSet(self, m, g_or_u, lmax); CHKERRQ(ierr);

  } else {

    int L; 
    PetscOptionsGetInt(NULL, NULL, "-y1s_L", &L, &find);
    if(!find)
      SETERRQ(self->comm, 1, "-y1s_lmax or -y1s_L is necessary");

    ierr = Y1sSetOne(self, L, m); CHKERRQ(ierr);
  }
  
  return 0;
}

PetscErrorCode Y1sGetSize(Y1s self, int *n) {
  *n = self->num;
  return 0;
}
PetscErrorCode Y1sGetMaxL(Y1s self, int *lmax) {
  *lmax = 0;
  for(int i = 0; i < self->num; i++) 
    if(self->ls[i] > *lmax)
      *lmax = self->ls[i];
  return 0;
}
PetscErrorCode Y1sGetLM(Y1s self, int i, int *l, int *m) {

  PetscErrorCode ierr;
  int n;
  ierr =Y1sGetSize(self, &n); CHKERRQ(ierr);
  if(i < 0 || n <= i) {
    SETERRQ(self->comm, 1, "index out of range");
  }

  *l = self->ls[i];
  *m = self->m;
  
  return 0;
}

PetscErrorCode Y1sCreateY1Mat(Y1s self, Mat *M) {
  
  PetscErrorCode ierr;
  ierr = Y1sCreateY1MatOther(self, self, M); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode Y1sCreateY1MatOther(Y1s self, Y1s other, Mat *M) {
    
  PetscErrorCode ierr;
  int na = self->num; 
  int nb = other->num; 
  ierr = MatCreate(self->comm, M);
  CHKERRQ(ierr);
  ierr = MatSetSizes(*M, PETSC_DECIDE, PETSC_DECIDE, na, nb);
  CHKERRQ(ierr);
  ierr = MatSetUp(*M);
  CHKERRQ(ierr);
  return 0;
  
}
PetscErrorCode Y1sSY1Mat(Y1s self, Mat M) {

  PetscLogEvent EVENT_id;
  PetscLogEventRegister("Y1sSY1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);  
  
  int n = self->num;
  for (int i = 0; i < n; i++) {
      MatSetValue(M, i, i, 1.0, INSERT_VALUES);
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);


  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;

}
PetscErrorCode Y1sLambdaY1Mat(Y1s self, Mat M) {

  PetscLogEvent EVENT_id;
  PetscLogEventRegister("Y1sLambdaY1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);  
  
  int n = self->num;
  for (int i = 0; i < n; i++) {
    int L = self->ls[i];
    if(L != 0)
      MatSetValue(M, i, i, 1.0*L*(L+1), INSERT_VALUES);
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;

}
PetscErrorCode Y1sPqY1Mat(Y1s self, int q, Mat M, PetscBool *non0) {

  PetscLogEvent EVENT_id;
  PetscLogEventRegister("Y1sPqY1Mat", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);

  int n = self->num;
  int *ls = self->ls;
  int m = self->m;
  PetscBool find = PETSC_FALSE;
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      PetscReal v = Y1ElePq(ls[i], q, ls[j], 
			    m,        m);
      if(fabs(v) > FLT_EPSILON) {
	find = PETSC_TRUE;
	MatSetValue(M, i, j, v, INSERT_VALUES);
      }
    }
  }

  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}
PetscErrorCode Y1sYqkY1MatOther(Y1s self, Y1s other, int q, int k, Mat M, PetscBool *non0) {

  PetscLogEvent EVENT_id;
  PetscLogEventRegister("Y1sYqkY1MatOther", 0, &EVENT_id);
  PetscLogEventBegin(EVENT_id, 0,0,0,0);
  

  PetscErrorCode ierr;

  int na = self->num;
  int *lsa = self->ls;
  int ma = self->m;

  int nb = other->num;
  int *lsb = other->ls;
  int mb = other->m;  
  
  PetscBool find = PETSC_FALSE;
  
  for(int i = 0; i < na; i++) {
    for(int j = 0; j < nb; j++) {
      PetscReal v = Y1EleYqk(lsa[i], q, lsb[j], 
			     ma,     k, mb);
      if(fabs(v) > FLT_EPSILON) {
	find = PETSC_TRUE;
	ierr = MatSetValue(M, i, j, v, INSERT_VALUES);
	CHKERRQ(ierr);
      }
    }
  }

  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }

  PetscLogEventEnd(EVENT_id, 0,0,0,0);
  return 0;
}

