#include <float.h>
#include <rescol/y2s.h>
#include <rescol/y1s.h>


// ---- Y2 Functions ----
Y2 Y2Exchange(Y2 a) {
  Y2 b = {a.l2, a.l1, a.l, a.m};
  return b;
}
PetscReal Y2RedY1q(Y2 yp, int q, Y2 y) {
  if(yp.l2 != y.l2)
    return 0.0;
  PetscReal t1 = sign(y.l + yp.l1 + y.l2 + q);
  PetscReal t2 = sqrt(1.0 * (2*y.l+1) * (2*yp.l+1));
  PetscReal t3 = Y1RedYq(yp.l1, q, y.l1);
  PetscReal t4 = wigner6j(yp.l1, yp.l, y.l2,
			  y.l,   y.l1, q);
  return t1*t2*t3*t4;
}
PetscReal Y2EleY1q(Y2 yp, int q, Y2 y) {
  PetscReal t1 = sign(yp.l - yp.m);
  PetscReal t2 = Y2RedY1q(yp, q, y);
  PetscReal t3 = wigner3j(yp.l,  q, y.l,
			  -yp.m, 0, y.m);
  return t1*t2*t3;
}
PetscReal Y2EleYYq(Y2 yp, int q, Y2 y) {

  if(yp.l != y.l)
    return 0;
  if(yp.m != y.m)
    return 0;

  PetscReal t1 = sign(yp.l2 + y.l1 + y.l);
  PetscReal t2 = Y1RedYq(yp.l2, q, y.l2);
  PetscReal t3 = Y1RedYq(yp.l1, q, y.l1);
  PetscReal t4 = wigner6j(yp.l1, yp.l2, yp.l,
			  y.l2,  y.l1,  q);

  return t1*t2*t3*t4;
}
PetscReal Y2ElePq1A(Y2 yp, int q, Y2 y) {
  return sqrt(4.0*M_PI/(2*q+1))*Y2EleY1q(yp, q, y);
}
PetscReal Y2ElePq2A(Y2 yp, int q, Y2 y) {
  Y2 yyp = Y2Exchange(yp);
  Y2 yy = Y2Exchange(y);
  return Y2ElePq1A(yyp, q, yy);
}
PetscReal Y2ElePq12(Y2 yp, int q, Y2 y) {
  return 4.0*M_PI/(2*q+1)*Y2EleYYq(yp, q, y);
}


// ---- Y2s Methods -----
PetscErrorCode Y2sCreate(MPI_Comm comm, Y2s *p_self) {
  Y2s self;
  PetscNew(&self);
  self->comm = comm;
  *p_self = self;
  return 0;
}
PetscErrorCode Y2sDestroy(Y2s *p_self) {

  PetscFree((*p_self)->y2_list);
  PetscFree(*p_self);
  return 0;

}

PetscErrorCode Y2sView(Y2s self, PetscViewer v) {

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

  if(iascii) {
    PetscViewerASCIIPrintf(v, "Y2s object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "num: %d\n", self->num);    
    PetscViewerASCIIPrintf(v, "Y[(L1, L2), L, M] = {\n");
    PetscViewerASCIIPushTab(v);
    for(int i = 0; i < self->num; i++) {
      Y2 y = self->y2_list[i];
      PetscViewerASCIIPrintf(v, "Y[(%d, %d), %d, %d],\n",
			     y.l1, y.l2, y.l, y.m);      
    }
    PetscViewerASCIIPopTab(v);
    PetscViewerASCIIPrintf(v, "}\n");    
    PetscViewerASCIIPopTab(v);
  } else if(isbinary) {

  } else if(isdraw) {

  }
  return 0;
}

PetscErrorCode Y2sSet(Y2s self, int m, int g_or_u, int p_or_m, int lmax) {

  if(m < 0)
    SETERRQ(self->comm, 1, "M must be 0 or positive");

  if(g_or_u != GERADE && g_or_u != UNGERADE)
    SETERRQ(self->comm, 1, "illegal g_or_u");

  int idx = 0;
  PetscMalloc1((lmax+1)*(lmax+1)*(lmax+1), &self->y2_list);
  for(int L = m; L <= lmax; L++)
    for(int L1 = m; L1 <= lmax; L1++) 
      for(int L2 = m; L2 <= lmax; L2++) {
	if( ((L%2==0) == (g_or_u == GERADE)) &&
	    ((L1+L2)%2==0) == (p_or_m == PLUS) &&
	    (L1+L2 <= lmax) &&
	    TriangleQ(L, L2, L1)) {
	  Y2 y2 = {L1, L2, L, m};
	  self->y2_list[idx] = y2;
	  idx++;
	}
      }
  self->num = idx;
  return 0;  
}
PetscErrorCode Y2sSetFromOptions(Y2s self) {

  char rot[10] = "sigma";
  char parity[10] = "gerade";
  char mirror[10] = "plus";
  int lmax = 2;
  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, "-y2s_rot", rot, 10, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-y2s_parity", parity, 10, NULL); 
  CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-y2s_mirror", mirror, 10, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-y2s_lmax", &lmax, NULL); CHKERRQ(ierr);
  
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

  int g_or_u = GERADE;
  if(strcmp(parity, "gerade") == 0)
    g_or_u = GERADE;
  else if(strcmp(parity, "ungerade") == 0)
    g_or_u = UNGERADE;
  else
    SETERRQ(self->comm, 1, "options -parity <- {gerade, ungerade}");

  int p_or_m = PLUS;
  if(strcmp(mirror, "plus") == 0)
    p_or_m = PLUS;
  else if(strcmp(mirror, "minus") == 0)
    p_or_m = MINUS;
  else
    SETERRQ(self->comm, 1, "options -mirror <- {plus, minus}");

  if(lmax < 0)
    SETERRQ(self->comm, 1, "options lmax must non negative integer");

  ierr = Y2sSet(self, m, g_or_u, p_or_m, lmax); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode Y2sGetSize(Y2s self, int *n) {
  
  *n = self->num;
  return 0;

}
PetscErrorCode Y2sGetMaxL(Y2s self, int *L) {
  int acc = 0;
  int n; Y2sGetSize(self, &n);
  for(int i = 0; i < n; i++) {
    if(acc < self->y2_list[i].l1)
      acc = self->y2_list[i].l1;
  }
  *L = acc;
  return 0;
}

PetscErrorCode Y2sCreateY2Mat(Y2s self, Mat *M) {

  int n; Y2sGetSize(self, &n);
  MatCreate(self->comm, M);
  MatSetSizes(*M, n, n, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode Y2sCreateY2Vec(Y2s self, Vec *v) {
  int n; Y2sGetSize(self, &n);
  VecCreate(self->comm, v);
  VecSetSizes(*v, PETSC_DECIDE, n);
  VecSetFromOptions(*v);
  return 0;
}
PetscErrorCode Y2sSY2Mat(Y2s self, Mat M) {

   int n; Y2sGetSize(self, &n);
   for(int i = 0; i < n; i++) {
     MatSetValue(M, i, i, 1.0, INSERT_VALUES);
   }
   MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
   return 0;

}
PetscErrorCode Y2sGuessY2Vec(Y2s self, int L1, int L2, Vec v) {

  int n; Y2sGetSize(self, &n);
  for(int i = 0; i < n; i++) 
    if(self->y2_list[i].l1 == L1 && 
       self->y2_list[i].l2 == L2) {
      VecSetValue(v, i, 1.0, INSERT_VALUES);
    }

  VecAssemblyBegin(v);
  VecAssemblyEnd(v);  
  return 0;

}
PetscErrorCode Y2sLambda1Y2Mat(Y2s self, Mat M, PetscBool *non0) {

   int n; Y2sGetSize(self, &n);

   PetscBool find = PETSC_FALSE;
   for(int i = 0; i < n; i++) {
     int l1 = self->y2_list[i].l1;
     if(l1 != 0) {
       find = PETSC_TRUE;
       PetscScalar v = 1.0*l1*(l1+1);
       MatSetValue(M, i, i, v, INSERT_VALUES);
     }
   }

  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }

  return 0;

}
PetscErrorCode Y2sLambda2Y2Mat(Y2s self, Mat M, PetscBool *non0) {

  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) {
    int l2 = self->y2_list[i].l2;
    if(l2 != 0) {
      find = PETSC_TRUE;
      PetscScalar v = 1.0*l2*(l2+1);
      MatSetValue(M, i, i, v, INSERT_VALUES);
    }
  }
  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }
  
  return 0;
}
PetscErrorCode Y2sPq1AY2Mat(Y2s self, int q, Mat M, PetscBool *non0) {
  
  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq1A(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > FLT_EPSILON) {
	find = PETSC_TRUE;
	MatSetValue(M, i, j, v, INSERT_VALUES);
      }
    }

  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }
  return 0;

}
PetscErrorCode Y2sPq2AY2Mat(Y2s self, int q, Mat M, PetscBool *non0) {

  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq2A(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > FLT_EPSILON * 10) {
	find = PETSC_TRUE;
	MatSetValue(M, i, j, v, INSERT_VALUES);
      }
    }
  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }
  return 0;

}
PetscErrorCode Y2sPq12Y2Mat(Y2s self, int q, Mat M, PetscBool *non0) {

  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq12(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > 4.0*FLT_EPSILON) {
	find = PETSC_TRUE;
	MatSetValue(M, i, j, v, INSERT_VALUES);
      }
    }

  if(find) {
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    *non0 = PETSC_TRUE;
  } else {
    *non0 = PETSC_FALSE;
  }
  return 0;

}

