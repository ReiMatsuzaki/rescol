#include <float.h>
#include <rescol/angmoment.h>


// --- Utils ---
PetscReal wigner3j(int a, int b, int c, int d, int e, int f) {
  return gsl_sf_coupling_3j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}
PetscReal wigner6j(int a, int b, int c, int d, int e, int f) {
  return gsl_sf_coupling_6j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}
PetscBool TriangleQ(int j1, int j2, int j3) {
  return (
	  (abs(j1-j2) <= j3 && j3 <= j1+j2) ||
	  (abs(j2-j3) <= j1 && j1 <= j2+j3) ||
	  (abs(j3-j1) <= j2 && j2 <= j3+j1)
	  );
  return 0;
}
int sign(int n) {
  if(n%2==0)
    return 1;
  else
    return -1;
}

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
PetscErrorCode Y1sCreate(Y1s *y1s, MPI_Comm comm) {

  PetscErrorCode ierr;
  Y1s _y1s;
  ierr = PetscNew(&_y1s); CHKERRQ(ierr);

  _y1s->comm = comm;
  *y1s = _y1s;
  return 0;
}
PetscErrorCode Y1sDestroy(Y1s* y1s) {

  PetscFree((*y1s)->ls);
  PetscFree(*y1s);
  return 0;

}
PetscErrorCode Y1sSet(Y1s y1s, int m, int g_or_u, int lmax) {

  PetscErrorCode ierr;
  if(m < 0)
    SETERRQ(y1s->comm, 1, "M must be 0 or positive");

  if(g_or_u != GERADE && g_or_u != UNGERADE)
    SETERRQ(y1s->comm, 1, "illegal g_or_u");
	    
  y1s->num = (lmax-m)/2 + 1;
  ierr = PetscMalloc1(y1s->num, &y1s->ls); CHKERRQ(ierr);
  y1s->m = m;

  int i = 0;
  for(int L = m; L <= lmax; L+=2) {
    y1s->ls[i] = L; i++;
  }

  return 0;
}
PetscErrorCode Y1sCreateFromOptions(Y1s *y1s, MPI_Comm comm) {

  char rot[10] = "sigma";
  char parity[10] = "gerade";
  int lmax = 2;
  PetscErrorCode ierr;
  ierr = PetscOptionsGetString(NULL, "-y1s_rot", rot, 10, NULL); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-y1s_parity", parity, 10, NULL); 
  CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, "-y1s_lmax", &lmax, NULL); CHKERRQ(ierr);
  
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
    SETERRQ(comm, 1, "options -rot <- {sigma, pi, delta, phi}");

  int g_or_u = GERADE;
  if(strcmp(parity, "gerade") == 0)
    g_or_u = GERADE;
  else if(strcmp(parity, "ungerade") == 0)
    g_or_u = UNGERADE;
  else
    SETERRQ(comm, 1, "options -parity <- {gerade, ungerade}");

  if(lmax < 0)
    SETERRQ(comm, 1, "options lmax must non negative integer");

  ierr = Y1sCreate(y1s, comm); CHKERRQ(ierr);
  ierr = Y1sSet(*y1s, m, g_or_u, lmax); CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode Y1sView(Y1s ys1) {

  PetscPrintf(ys1->comm, "num: %d\n", ys1->num);
  PetscPrintf(ys1->comm, "m: %d\n", ys1->m);
  PetscPrintf(ys1->comm, "ls: ");
  for(int i = 0; i < ys1->num; i++) {
    PetscPrintf(ys1->comm, "%d ", ys1->ls[i]);
  }
  PetscPrintf(ys1->comm, "\n");
  return 0;

}
PetscErrorCode Y1sGetSize(Y1s y1s, int *n) {
  *n = y1s->num;
  return 0;
}
PetscErrorCode Y1sGetMaxL(Y1s self, int *lmax) {
  *lmax = 0;
  for(int i = 0; i < self->num; i++) 
    if(self->ls[i] > *lmax)
      *lmax = self->ls[i];
  return 0;
}
PetscErrorCode Y1sCreateY1Mat(Y1s self, Mat *M) {
  
  int n = self->num;
  MatCreate(self->comm, M);
  MatSetSizes(*M, n, n, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode Y1sSetSY1Mat(Y1s self, Mat *M) {
  
  Y1sCreateY1Mat(self, M);
  int n = self->num;
  for (int i = 0; i < n; i++) {
      MatSetValue(*M, i, i, 1.0, INSERT_VALUES);
  }
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;

}
PetscErrorCode Y1sSetLambdaY1Mat(Y1s y1s, Mat *M) {

  Y1sCreateY1Mat(y1s, M);
  int n = y1s->num;
  for (int i = 0; i < n; i++) {
    int L = y1s->ls[i];
    if(L != 0)
      MatSetValue(*M, i, i, 1.0*L*(L+1), INSERT_VALUES);
  }
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;

}
PetscErrorCode Y1sSetPqY1Mat(Y1s y1s, int q, Mat *M) {

  Y1sCreateY1Mat(y1s, M);
  int n = y1s->num;
  int *ls = y1s->ls;
  int m = y1s->m;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      PetscReal v = Y1ElePq(ls[i], q, ls[j], 
			    m,        m);
      MatSetValue(*M, i, j, v, INSERT_VALUES);
    }
  }
  MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  return 0;
}

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
PetscErrorCode Y2sCreate(Y2s *y2s, MPI_Comm comm) {
  Y2s _y2s;
  PetscNew(&_y2s);
  _y2s->comm = comm;
  *y2s = _y2s;
  return 0;
}
PetscErrorCode Y2sDestroy(Y2s *y2s) {

  PetscFree((*y2s)->y2_list);
  PetscFree(*y2s);
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
PetscErrorCode Y2sCreateFromOptions(Y2s *y2s, MPI_Comm comm) {
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
    SETERRQ(comm, 1, "options -rot <- {sigma, pi, delta, phi}");

  int g_or_u = GERADE;
  if(strcmp(parity, "gerade") == 0)
    g_or_u = GERADE;
  else if(strcmp(parity, "ungerade") == 0)
    g_or_u = UNGERADE;
  else
    SETERRQ(comm, 1, "options -parity <- {gerade, ungerade}");

  int p_or_m = PLUS;
  if(strcmp(mirror, "plus") == 0)
    p_or_m = PLUS;
  else if(strcmp(mirror, "minus") == 0)
    p_or_m = MINUS;
  else
    SETERRQ(comm, 1, "options -mirror <- {plus, minus}");

  if(lmax < 0)
    SETERRQ(comm, 1, "options lmax must non negative integer");

  ierr = Y2sCreate(y2s, comm); CHKERRQ(ierr);
  ierr = Y2sSet(*y2s, m, g_or_u, p_or_m, lmax); CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode Y2sView(Y2s self) {

  int n; Y2sGetSize(self, &n);
  PetscPrintf(self->comm, "num: %d\n", n);
  PetscPrintf(self->comm, "(L1, L2, L, M)\n");
  for(int i = 0; i < n; i++) {
    Y2 y = self->y2_list[i];
    PetscPrintf(self->comm, "(%d, %d, %d, %d)\n",
		y.l1, y.l2, y.l, y.m);
  }
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
PetscErrorCode Y2sInitY2Mat(Y2s self, Mat *M) {

  int n; Y2sGetSize(self, &n);
  MatCreate(self->comm, M);
  MatSetSizes(*M, n, n, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;

}
PetscErrorCode Y2sSetSY2Mat(Y2s self, Mat *M) {

   Y2sInitY2Mat(self, M);
   int n; Y2sGetSize(self, &n);
   for(int i = 0; i < n; i++) {
     MatSetValue(*M, i, i, 1.0, INSERT_VALUES);
   }
   MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
   return 0;

}
PetscErrorCode Y2sSetGuessY2Vec(Y2s self, int L1, int L2, Vec *v) {

  VecCreate(self->comm, v);
  int n; Y2sGetSize(self, &n);
  VecSetSizes(*v, PETSC_DECIDE, n);
  VecSetFromOptions(*v);
  
  for(int i = 0; i < n; i++) 
    if(self->y2_list[i].l1 == L1 && 
       self->y2_list[i].l2 == L2) {
      VecSetValue(*v, i, 1.0, INSERT_VALUES);
    }

  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);  
  return 0;

}
PetscErrorCode Y2sSetLambda1Y2Mat(Y2s self, Mat *M) {

   Y2sInitY2Mat(self, M);
   int n; Y2sGetSize(self, &n);

   PetscBool find = PETSC_FALSE;
   for(int i = 0; i < n; i++) {
     int l1 = self->y2_list[i].l1;
     if(l1 != 0) {
       find = PETSC_TRUE;
       PetscScalar v = 1.0*l1*(l1+1);
       MatSetValue(*M, i, i, v, INSERT_VALUES);
     }
   }

   if(find) {
     MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
   } else {
     MatDestroy(M);
     *M = NULL;
   }
   return 0;

}
PetscErrorCode Y2sSetLambda2Y2Mat(Y2s self, Mat *M) {

  Y2sInitY2Mat(self, M);
  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) {
    int l2 = self->y2_list[i].l2;
    if(l2 != 0) {
      find = PETSC_TRUE;
      PetscScalar v = 1.0*l2*(l2+1);
      MatSetValue(*M, i, i, v, INSERT_VALUES);
    }
  }
  if(find) {
    MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  } else {
    MatDestroy(M);
    *M = NULL;
  }
  
  return 0;
}
PetscErrorCode Y2sSetPq1AY2Mat(Y2s self, int q, Mat *M) {

  Y2sInitY2Mat(self, M);
  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq1A(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > FLT_EPSILON) {
	find = PETSC_TRUE;
	MatSetValue(*M, i, j, v, INSERT_VALUES);
      }
    }

  if(find) {
    MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  } else {
    MatDestroy(M);
    *M = NULL;
  }
  return 0;

}
PetscErrorCode Y2sSetPq2AY2Mat(Y2s self, int q, Mat *M) {

  Y2sInitY2Mat(self, M);
  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq2A(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > FLT_EPSILON * 10) {
	find = PETSC_TRUE;
	MatSetValue(*M, i, j, v, INSERT_VALUES);
      }
    }
  if(find) {
    MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  } else {
    MatDestroy(M);
    *M = NULL;
  }
  return 0;

}
PetscErrorCode Y2sSetPq12Y2Mat(Y2s self, int q, Mat *M) {

  Y2sInitY2Mat(self, M);
  int n; Y2sGetSize(self, &n);

  PetscBool find = PETSC_FALSE;
  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++) {
      PetscReal v = Y2ElePq12(self->y2_list[i], q, self->y2_list[j]);
      if(fabs(v) > FLT_EPSILON) {
	find = PETSC_TRUE;
	MatSetValue(*M, i, j, v, INSERT_VALUES);
      }
    }

  if(find) {
    MatAssemblyBegin(*M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*M, MAT_FINAL_ASSEMBLY);
  } else {
    MatDestroy(M);
    *M = NULL;
  }
  return 0;

}