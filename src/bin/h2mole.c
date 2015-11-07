#include <rescol/oce2.h>

#if defined(PETSC_USE_COMPLEX)
#define ABS cabs
#else
#define ABS fabs
#endif

static char help[] = "solve H2 mole eigen value problem";
/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  -eri
  -bond_length : 

  -y2s_rot
  -y2s_parity
  -y2s_mirror
  -y2s_lmax

  -fem_type
  -fd_num  -fd_xmax
  -bss_order 
  -dvr_nq
  -bps_type -bps_num_zs -bps_zmax

  -in_dir -out_dir :
*/

struct _p_H2 {

  // -----  Constant -----
  MPI_Comm comm;
  char in_dir[100];
  char out_dir[100];
  char guess_type[10];
  char eri_option[10];

  OCE2 oce2;

  Mat H0;
  Mat S;
  
  PetscReal d_bondlength; // delta of bond length
  PetscInt num_bondlength; // number of bond length for caluclation

  // ---- varable -----
  PetscReal bondlength;  // initial bond length
  Mat H;
  Vec *guess;
  int num_guess; // # of guess vector to use
  int max_num_guess; // size of array "guess"
};

typedef struct _p_H2* H2;

PetscErrorCode H2CreateFromOptions(H2 *h2, MPI_Comm comm) {

  PrintTimeStamp(comm, "Init", NULL);
  PetscErrorCode ierr;
  H2 _h2;
  ierr = PetscMalloc1(1, &_h2); CHKERRQ(ierr);
  *h2 = NULL;

  _h2->comm = comm;
  strcpy(_h2->in_dir, ".");
  strcpy(_h2->out_dir, ".");
  strcpy(_h2->guess_type, ".");
  strcpy(_h2->eri_option, ".");

  _h2->bondlength = 2.0;
  _h2->d_bondlength = 1.0;
  _h2->num_bondlength = 1;

  PetscOptionsGetString(NULL, "-in_dir", _h2->in_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-out_dir", _h2->out_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-eri", _h2->eri_option, 10, NULL);  
  PetscOptionsGetString(NULL, "-guess_type", _h2->guess_type, 10, NULL);
  PetscOptionsGetReal(NULL, "-bondlength", &_h2->bondlength, NULL);
  PetscOptionsGetReal(NULL, "-d_bondlength", &_h2->d_bondlength, NULL);
  PetscOptionsGetInt(NULL, "-num_bondlength", &_h2->num_bondlength, NULL);

  ierr = OCE2CreateFromOptions(&_h2->oce2, comm); CHKERRQ(ierr);

  if(_h2->d_bondlength <= 0.0)
    SETERRQ(comm, 1, "d_bondlength > 0");
  if(_h2->num_bondlength < 1)
    SETERRQ(comm, 1, "num_bondlength > 0");

  _h2->H0 = NULL;
  _h2->S = NULL;

  PetscOptionsGetReal(NULL, "-bondlength0", &_h2->bondlength, NULL);
  _h2->H = NULL;
  _h2->num_guess = 1;
  int nev = 0; PetscOptionsGetInt(NULL, "-eps_nev", &nev, NULL);
  if(nev > 0)
    _h2->max_num_guess = nev;
  else
    _h2->max_num_guess = 1;
  ierr = PetscMalloc1(_h2->max_num_guess, &_h2->guess);
  for(int i = 0; i < _h2->max_num_guess; i++)
    _h2->guess[i] = NULL;

  *h2 = _h2;
  return 0;
}
PetscErrorCode H2Destroy(H2 *h2) {
  PetscErrorCode ierr;
  H2 this = *h2;

  ierr = OCE2Destroy(&this->oce2);

  ierr = MatDestroy(&this->H0); CHKERRQ(ierr);
  if(this->S != NULL) {
    ierr = MatDestroy(&this->S); CHKERRQ(ierr);
  }
  
  ierr = MatDestroy(&this->H); CHKERRQ(ierr);

  for(int i = 0; i < this->max_num_guess; i++) {
    ierr = VecDestroy(&this->guess[i]); CHKERRQ(ierr);
  }
  ierr = PetscFree(this->guess); CHKERRQ(ierr);

  ierr = PetscFree(*h2); CHKERRQ(ierr);
  
  return 0;
}
PetscErrorCode H2SetBasic(H2 this) {

  PrintTimeStamp(this->comm, "Mat0", NULL);
  PetscErrorCode ierr;
 
   // guess
   if(strcmp(this->guess_type, "read") == 0) {
     PrintTimeStamp(this->comm, "read guess", NULL);
     VecCreate(this->comm, &this->guess[0]);

     PetscViewer viewer;
     char path[100]; sprintf(path, "%s/guess.vec.dat", this->in_dir);
     ierr = PetscViewerBinaryOpen(this->comm, path, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
     ierr = VecLoad(this->guess[0], viewer); CHKERRQ(ierr);        
     ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

   } else
     SETERRQ(this->comm, 1, "guess_type<-{read}");

   return 0;
 }
PetscErrorCode H2SetH0_and_S(H2 this) {
  PetscErrorCode ierr;

  PrintTimeStamp(this->comm, "TMat", NULL);
  ierr = OCE2SetTMat(this->oce2, &this->H0); CHKERRQ(ierr);

  
  if(strcmp(this->eri_option, "direct") == 0) {
    PrintTimeStamp(this->comm, "Vee", NULL);
    ierr = OCE2PlusVeeMat(this->oce2, &this->H0); CHKERRQ(ierr);
  }
  
  PetscBool s_is_id; FEMInfGetOverlapIsId(this->oce2->fem, &s_is_id);
  if(!s_is_id) {
    PrintTimeStamp(this->comm, "S", NULL);
    ierr = OCE2SetSMat(this->oce2, &this->S); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode H2SetH_as_H0_plus_NE(H2 this) {

  PetscErrorCode ierr;

  PrintTimeStamp(this->comm, "NE", NULL);
  PetscPrintf(this->comm, "bond_length: %f\n", this->bondlength);
  
  PetscScalar a = this->bondlength/2.0;
  if(this->H == NULL) {
    ierr = MatConvert(this->H0, MATSAME, MAT_INITIAL_MATRIX, &this->H);
  } else {
    ierr = MatCopy(this->H0, this->H, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  }

  ierr = OCE2PlusVneMat(this->oce2, a, 1.0, &this->H); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode CheckSymmetry(H2 this, Vec cs, PetscBool *is_sym) {
  // assuming size of cs is equal to (n_r1)(n_r1)(n_y2)

  PetscInt n_r1, n_y2, size_cs;
  OCE2GetSizes(this->oce2, &n_r1, &n_y2);
  VecGetSize(cs, &size_cs);

  PetscPrintf(this->comm, "Vec size: %d\n", size_cs);
  PetscPrintf(this->comm, "r1: %d\n", n_r1);
  PetscPrintf(this->comm, "y2: %d\n", n_y2);
  
  PetscScalar *vs;
  PetscReal eps = 0.000001;
  PetscMalloc1(2, &vs);
  for(int y = 0; y < n_y2; y++)
    for(int i = 0; i < n_r1; i++) {
      for(int j = 0; j < i; j++) {
	PetscInt idx[2]; 
	idx[0] = i + j*n_r1 + y*n_r1*n_r1; 
	idx[1] = j + i*n_r1 + y*n_r1*n_r1;
	VecGetValues(cs, 2, idx, vs);
	if(ABS(vs[0]) > eps && ABS(vs[1]) > eps)
	  goto end;
    }
  }

 end:
  *is_sym = (ABS(vs[0]*vs[1]) > 0.0);
  if(fabs(ABS(vs[0])-ABS(vs[1])) > eps | ABS(vs[0]) < eps) {
    PetscPrintf(this->comm, "WARNING: check symmetry may contain error\n");
    PetscPrintf(this->comm, "vs[0] = %f\n", vs[0]);
    PetscPrintf(this->comm, "vs[1] = %f\n", vs[1]);
  }

  PetscFree(vs);
  return 0;
}
PetscErrorCode H2Solve(H2 this) {

  PrintTimeStamp(this->comm, "EPS", NULL);
  MPI_Comm comm = this->comm;
  PetscErrorCode ierr;
  EPS eps;

  ierr = EPSCreate(this->comm, &eps);  CHKERRQ(ierr);
  ierr = EPSSetTarget(eps, -4.0);  CHKERRQ(ierr);

  PetscBool s_is_id; FEMInfGetOverlapIsId(this->oce2->fem, &s_is_id);
  if(s_is_id) 
    EPSSetProblemType(eps, EPS_HEP); 
  else
    EPSSetProblemType(eps, EPS_GHEP); 

  if(this->guess[0] != NULL) {
    ierr = EPSSetInitialSpace(eps, this->num_guess, this->guess); CHKERRQ(ierr);
  }
  
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 
  ierr = EPSSetOperators(eps, this->H, this->S); CHKERRQ(ierr);
  EPSSetFromOptions(eps);
  PrintTimeStamp(this->comm, "solve", NULL);
  
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  PetscPrintf(comm, "==== Results ====\n");
  int nconv;
  PetscScalar kr;
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    if(this->guess[i] == NULL) {
      ierr = MatCreateVecs(this->H, NULL, &this->guess[i]); CHKERRQ(ierr);
    }
    EPSGetEigenpair(eps, i, &kr, NULL, this->guess[i], NULL);
    PetscBool is_sym;
    CheckSymmetry(this, this->guess[i], &is_sym);
    PetscPrintf(comm, "eig%i: %f, (%s)\n", i, kr, is_sym?"sym":"anti");
  }
  this->num_guess = nconv;
  EPSDestroy(&eps);
  return 0;
}
PetscErrorCode H2DestroyMini(H2 this) {
  PrintTimeStamp(this->comm, "Destroy", NULL);
  MatDestroy(&this->H);
  return 0;
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  H2 h2;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  PetscPrintf(comm, "Scalar: complex\n");
#else
  PetscPrintf(comm, "Scalar: real\n");
#endif
  PetscOptionsBegin(comm, "", "h2mole.c options", "none");
  ierr = H2CreateFromOptions(&h2, comm); CHKERRQ(ierr);
  PetscOptionsEnd();

  OCE2View(h2->oce2);
  
  ierr = H2SetBasic(h2); CHKERRQ(ierr);
  ierr = H2SetH0_and_S(h2); CHKERRQ(ierr);

  for(int i = 0; i < h2->num_bondlength; i++) {
    PetscPrintf(comm, "\n");
    ierr = H2SetH_as_H0_plus_NE(h2); CHKERRQ(ierr);
    ierr = H2Solve(h2); CHKERRQ(ierr);
    ierr = H2DestroyMini(h2); CHKERRQ(ierr);

    h2->bondlength += h2->d_bondlength;
  }

  ierr = H2Destroy(&h2);CHKERRQ(ierr);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
