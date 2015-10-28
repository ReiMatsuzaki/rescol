#include "fem_inf.h"

static char help[] = "solve H2 mole eigen value problem";
/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  -bond_length : 
  -lmax :
  -gerade or -ungerade:
  -eri

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

  int qmax;
  FEMInf fem;
  Mat s_r1;
  Mat s_y2;
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

  _h2->qmax = 10;
  _h2->bondlength = 2.0;
  _h2->d_bondlength = 1.0;
  _h2->num_bondlength = 1;

  PetscOptionsGetString(NULL, "-in_dir", _h2->in_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-out_dir", _h2->out_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-guess_type", _h2->guess_type, 10, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &_h2->qmax, NULL);
  PetscOptionsGetString(NULL, "-eri", _h2->eri_option, 10, NULL);  
  PetscOptionsGetReal(NULL, "-d_bondlength", &_h2->d_bondlength, NULL);
  PetscOptionsGetInt(NULL, "-num_bondlength", &_h2->num_bondlength, NULL);

  ierr = FEMInfCreateFromOptions(&_h2->fem, comm); CHKERRQ(ierr);

  if(_h2->d_bondlength <= 0.0)
    SETERRQ(comm, 1, "d_bondlength > 0");
  if(_h2->num_bondlength < 1)
    SETERRQ(comm, 1, "num_bondlength > 0");

  _h2->s_r1 = NULL;
  _h2->s_y2 = NULL;
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

  ierr = FEMInfDestroy(&this->fem);

  ierr = MatDestroy(&this->s_r1); CHKERRQ(ierr);
  ierr = MatDestroy(&this->s_y2); CHKERRQ(ierr);

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
  PetscErrorCode ierr = FEMInfSetSR1Mat(this->fem, &this->s_r1); CHKERRQ(ierr); 
  ierr = MatSetDirFile(this->in_dir, "s_y2mat.dat", &this->s_y2); 
  CHKERRQ(ierr);

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

PetscErrorCode H2SetSMat(H2 this) {

  PrintTimeStamp(this->comm, "SMat", NULL);
  PetscErrorCode ierr;
  ierr = MatSetSynthesize3Fast(this->s_r1, this->s_r1, this->s_y2, 
			       this->comm, &this->S);
  return 0;
}

PetscErrorCode H2SetD2(H2 this) {

  PrintTimeStamp(this->comm, "D2Mat", NULL);
  PrintTimeStamp(this->comm, "D2(r1)", NULL);
  
  // D2 
  Mat d2_r1;
  PetscErrorCode ierr;
  ierr = FEMInfSetD2R1Mat(this->fem, &d2_r1); CHKERRQ(ierr);
  ierr = MatScale(d2_r1, -0.5); CHKERRQ(ierr);

  Mat D;
  PrintTimeStamp(this->comm, "D2(syn)", NULL);
  ierr = MatSetSynthesize3Fast(this->s_r1, d2_r1, this->s_y2, 
			       this->comm, &this->H0); 
  CHKERRQ(ierr);
  ierr = MatSetSynthesize3Fast(d2_r1, this->s_r1, this->s_y2, 
			       this->comm, &D); CHKERRQ(ierr);

  PrintTimeStamp(this->comm, "D2(axpy)", NULL);
  ierr = MatAXPY(this->H0, 1.0, D, DIFFERENT_NONZERO_PATTERN);  CHKERRQ(ierr);

  PrintTimeStamp(this->comm, "D2(free)", NULL);
  MatDestroy(&d2_r1); MatDestroy(&D);
  return 0;
}

PetscErrorCode H2AddL(H2 this) {

  PrintTimeStamp(this->comm, "LMat", NULL);

  Mat l_r1; 
  PetscErrorCode ierr;
  ierr = FEMInfSetR2invR1Mat(this->fem, &l_r1); CHKERRQ(ierr);
  
  char l1_path[100]; sprintf(l1_path, "%s/l_1_y2mat.dat", this->in_dir); 
  FILE *fl1 = fopen(l1_path, "r");
  if(fl1 != NULL) {
    Mat l_1_y2, L;
    ierr = MatCreateFromCOOFormatFileHandler(fl1, &l_1_y2); CHKERRQ(ierr);
    ierr = MatSetSynthesize3Fast(l_r1, this->s_r1, l_1_y2, this->comm, &L); 
    CHKERRQ(ierr);
    MatAXPY(this->H0, 0.5, L, SUBSET_NONZERO_PATTERN);
    MatDestroy(&L); MatDestroy(&l_1_y2); fclose(fl1);
  }

  char l2_path[100]; sprintf(l2_path, "%s/l_2_y2mat.dat", this->in_dir); 
  FILE *fl2 = fopen(l2_path, "r");
  if(fl2 != NULL) {
    Mat l_2_y2, L;
    ierr = MatCreateFromCOOFormatFileHandler(fl1, &l_2_y2); CHKERRQ(ierr);
    ierr = MatSetSynthesize3Fast(this->s_r1, l_r1, l_2_y2, this->comm, &L);
    CHKERRQ(ierr);
    MatAXPY(this->H0, 0.5, L, SUBSET_NONZERO_PATTERN);
    MatDestroy(&L); MatDestroy(&l_2_y2); fclose(fl2);
  }

  MatDestroy(&l_r1);
  return 0;
}

PetscErrorCode H2AddVee(H2 this) {

  PrintTimeStamp(this->comm, "EE", NULL);

  PetscErrorCode ierr;
  for(int q = 0; q < this->qmax; q++) {
    char path[100]; sprintf(path, "%s/p%d_12_y2mat.dat", this->in_dir, q);
    FILE *f = fopen(path, "r"); 
    if(f != NULL) {
      PetscPrintf(this->comm, "%d ", q);
      Mat r2, y2, V; 
      ierr = FEMInfSetEER2Mat(this->fem, q, &r2); CHKERRQ(ierr);
      ierr = MatCreateFromCOOFormatFileHandler(f, &y2); CHKERRQ(ierr);
      ierr = MatSetSynthesizeFast(r2, y2, this->comm, &V); CHKERRQ(ierr);
      ierr = MatAXPY(this->H0, 1.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&r2); MatDestroy(&y2); MatDestroy(&V);
      fclose(f);
    }
  }
  PetscPrintf(this->comm, "\n");
  return 0;
}

PetscErrorCode H2SetH_as_H0_plus_NE(H2 this) {

  PetscErrorCode ierr;

  PrintTimeStamp(this->comm, "NE", NULL);
  PetscPrintf(this->comm, "bond_length: %f\n", this->bondlength);
  
  PetscScalar a = this->bondlength/2.0;
  if(this->H == NULL) {
    // ierr = MatDuplicate(this->H0, MAT_COPY_VALUES, &this->H); CHKERRQ(ierr);
    ierr = MatConvert(this->H0, MATSAME, MAT_INITIAL_MATRIX, &this->H);
  } else {
    ierr = MatCopy(this->H0, this->H, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  }

  for(int q = 0; q < this->qmax; q++) {
    char path1[100]; char path2[100]; 
    sprintf(path1, "%s/p%d_A1_y2mat.dat", this->in_dir, q);
    sprintf(path2, "%s/p%d_A2_y2mat.dat", this->in_dir, q);
    FILE* f1 = fopen(path1, "r"); FILE* f2 = fopen(path2, "r");
    if(f1 != NULL && f2 != NULL) {
      PetscPrintf(this->comm, "%d ", q);
      Mat q_r1;
      ierr = FEMInfSetENR1Mat(this->fem, q, a, &q_r1);  CHKERRQ(ierr);

      Mat pq_1_y2;
      ierr = MatCreateFromCOOFormatFileHandler(f1, &pq_1_y2); CHKERRQ(ierr);

      Mat V1;
      ierr = MatSetSynthesize3Fast(q_r1, this->s_r1, pq_1_y2, this->comm, &V1); CHKERRQ(ierr);
      MatAXPY(this->H, -2.0, V1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&pq_1_y2); MatDestroy(&V1);

      Mat pq_2_y2;
      ierr = MatCreateFromCOOFormatFileHandler(f2, &pq_2_y2); CHKERRQ(ierr);

      Mat V2;
      ierr = MatSetSynthesize3Fast(this->s_r1, q_r1, pq_2_y2, this->comm, &V2); 
      CHKERRQ(ierr);
      ierr = MatAXPY(this->H, -2.0, V2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      MatDestroy(&pq_2_y2); MatDestroy(&V2); 

      MatDestroy(&q_r1); fclose(f1); fclose(f2);
    }
  }
  PetscPrintf(this->comm, "\n");
  return 0;
}

PetscErrorCode H2Solve(H2 this) {

  PrintTimeStamp(this->comm, "EPS", NULL);
  MPI_Comm comm = this->comm;
  PetscErrorCode ierr;
  EPS eps;

  ierr = EPSCreate(this->comm, &eps);  CHKERRQ(ierr);
  ierr = EPSSetTarget(eps, -4.0);  CHKERRQ(ierr);

  PetscBool s_is_id; FEMInfGetOverlapIsId(this->fem, &s_is_id);
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
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
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
  PetscOptionsBegin(comm, "", "h2mole.c options", "none");
  ierr = H2CreateFromOptions(&h2, comm); CHKERRQ(ierr);
  PetscOptionsEnd();

  PetscBool s_is_id; FEMInfGetOverlapIsId(h2->fem, &s_is_id);
  FEMInfFPrintf(h2->fem, stdout, 0);
  
  ierr = H2SetBasic(h2); CHKERRQ(ierr);
  if(!s_is_id)
    ierr = H2SetSMat(h2); CHKERRQ(ierr);
  ierr = H2SetD2(h2); CHKERRQ(ierr);
  ierr = H2AddL(h2); CHKERRQ(ierr);
  if(strcmp(h2->eri_option, "direct") == 0) 
    ierr = H2AddVee(h2); CHKERRQ(ierr);

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
