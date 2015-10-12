#include "angmoment.h"

// --- Function ---
double PyGaunt(int j1, int j2, int j3, 
	       int m1, int m2, int m3) {

  PyObject *pName, *pModule, *pFunc, *pArgs, *pValue;
  Py_Initialize();
  //  pName = PyString_FromString("sympy.physics.wigner");
  pName = PyString_FromString("sympy");
  
  //  pModule = PyImport_ImportModule("sympy.physics.wigner");
  pModule = PyImport_Import(pName);
  if(pModule == NULL) {
    printf("failed to bind py module\n");
    exit(1);
  }
    //    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
    //	    

  pFunc = PyObject_GetAttrString(pModule, "gaunt");
  if(pFunc == NULL) 
    printf("failed to bind py func\n");
//    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
//	    "failed to bind py func");

  pArgs = PyTuple_New(7);
  pValue = PyInt_FromLong(j1); PyTuple_SetItem(pArgs, 0, pValue);
  pValue = PyInt_FromLong(j2); PyTuple_SetItem(pArgs, 1, pValue);
  pValue = PyInt_FromLong(j3); PyTuple_SetItem(pArgs, 2, pValue);
  pValue = PyInt_FromLong(m1); PyTuple_SetItem(pArgs, 3, pValue);
  pValue = PyInt_FromLong(m2); PyTuple_SetItem(pArgs, 4, pValue);
  pValue = PyInt_FromLong(m3); PyTuple_SetItem(pArgs, 5, pValue);
  pValue = PyInt_FromLong(10); PyTuple_SetItem(pArgs, 6, pValue);
  
  pValue = PyObject_CallObject(pFunc, pArgs);
  if(pValue == NULL) 
    printf("failed to bind py value");

  double res = PyFloat_AsDouble(pValue);
  Py_DECREF(pModule); Py_DECREF(pFunc); Py_DECREF(pArgs); Py_DECREF(pValue);

  Py_Finalize();
  return res;
}
double Y1RedMat_Yq(int j1, int j2, int j3) {
  return 1.0;
}
double Y1Mat_Yqk(int j1, int j2, int j3, int m1, int m2, int m3) {
  return 7.1;
}
// --- Y1s ---

PetscErrorCode Y1sCreate(Y1s *y1s, int l0, int maxl, Y1Sym sym, int m) {

  if(sym != GERADE && sym != UNGERADE) 
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, 
	    "sym must be GERADE or UNGERADE\n");     

  if(sym == GERADE && ((l0%2)==1 || (maxl%2)==1) )
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, 
	    "l0 and maxl must be even number under gerade symmetry");
  if(sym == UNGERADE && ((l0%2)==0 || (maxl%2)==0) )
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, 
	    "l0 and maxl must be odd number for UNGERADE symmetry");

  Y1s _y1s = malloc(sizeof(struct _p_Y1s));
  *y1s = NULL;

  _y1s->num = (maxl - l0)/2;
  _y1s->ls = (int*)malloc(sizeof(int)*_y1s->num);
  _y1s->m = m;

  int i = 0;
  for(int L = l0; L <= maxl; L+=2) {
    _y1s->ls[i] = L; i++;
  }

  *y1s = _y1s;
  return 0;
}

PetscErrorCode Y1sDestroy(Y1s* y1s) {
  free((*y1s)->ls);
  free(*y1s);
  return 0;
}

PetscErrorCode Y1sCreateY1Mat(Y1s this, MPI_Comm comm, Mat *M) {
  int n = this->num;
  MatCreate(comm, M);
  MatSetSizes(*M, n, n, n, n);
  MatSetFromOptions(*M);
  MatSetUp(*M);
  return 0;
}

PetscErrorCode Y1sCalcLambdaY1Mat(Y1s this, Mat M, InsertMode mode) {
  int n = this->num;
  for (int i = 0; i < n; i++) {
    int L = this->ls[i];
    if(L != 0)
      MatSetValue(M, i, i, 1.0*L*(L+1), mode);
  }
  return 0;
}

PetscErrorCode Y1sCalcPqY1Mat(Y1s this, int q, Mat M, InsertMode mode) {
  return 1;
}
