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

  *h2 = _h2;
  return 0;
}

PetscErrorCode H2Destroy(H2 *h2) {
  PetscErrorCode ierr;
  H2 this = *h2;

  ierr = FEMInfDestroy(&this->fem);

  if(this->S != NULL) {
    ierr = MatDestroy(&this->S); CHKERRQ(ierr);
  }
  ierr = MatDestroy(&this->H0); CHKERRQ(ierr);
  for(int i = 0; i < this->max_num_guess; i++) {
    ierr = VecDestroy(&this->guess[i]); CHKERRQ(ierr);
  }
  ierr = PetscFree(this->guess); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode H2SetBasic(H2 this) {

  PrintTimeStamp(this->comm, "Mat0", NULL);
  PetscErrorCode ierr = FEMInfSetSR1Mat(this->fem, &this->s_r1); CHKERRQ(ierr); 
  ierr = MatSetDirFile(this->in_dir, "s_y2mat.dat", &this->s_y2); 
  CHKERRQ(ierr);

  if(strcmp(this->guess_type, "read") == 0) {
    PrintTimeStamp(this->comm, "read guess", NULL);
    PetscViewer viewer;
    VecCreate(this->comm, &this->guess[0]);
    char path[100]; sprintf(path, "%s/guess.vec.dat", this->in_dir);
    PetscViewerBinaryOpen(this->comm, path, FILE_MODE_READ, &viewer);
    ierr = VecLoad(this->guess[0], viewer); CHKERRQ(ierr);        
  } else if (strcmp(this->guess_type, "none") == 0) {
    this->guess[0] = NULL;
  } else
    SETERRQ(this->comm, 1, "guess_type<-{read, none}");

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
  // ierr = MatConvert(this->H0, MATSAME, MAT_INITIAL_MATRIX, &H);
  ierr = MatDuplicate(this->H0, MAT_COPY_VALUES, &this->H);CHKERRQ(ierr);

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
  //  ierr = EPSSetOperators(eps, this->H, this->S); CHKERRQ(ierr);
  ierr = EPSSetTarget(eps, -4.0);  CHKERRQ(ierr);

  PetscBool s_is_id; FEMInfGetOverlapIsId(this->fem, &s_is_id);
  if(s_is_id) 
    EPSSetProblemType(eps, EPS_HEP); 
  else
    EPSSetProblemType(eps, EPS_GHEP); 

  if(this->guess[0] != NULL) {
    ierr = EPSSetInitialSpace(eps, 1, this->guess); CHKERRQ(ierr);
  }
  
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 

  ierr = EPSSetOperators(eps, this->H, this->S); CHKERRQ(ierr);
  PrintTimeStamp(this->comm, "solve", NULL);
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  PetscPrintf(comm, "==== Results ====\n");
  int nconv;
  PetscScalar kr;
  Vec xr;
  ierr = MatCreateVecs(this->H, NULL, &xr); CHKERRQ(ierr);
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, NULL, xr, NULL);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
    VecDuplicate(xr, &this->guess[i]);
  }  

  EPSDestroy(&eps);
  VecDestroy(&xr);
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

  return 0;
}
/*
int main(int argc, char **args) {

PetscErrorCode ierr;
  FEMInf fem;
  MPI_Comm comm = PETSC_COMM_SELF;
  char in_dir[100] = "default_in";
  char out_dir[100] = "default_out";
  char guess_type[10] = "none";
  double bond_length = 2.0;
  int qmax = 10;
  // none => not include ERI, direct=>directly evaluate ERI
  char eri_option[10] = "direct"; 

  // Initialize
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  time_t t0; PrintTimeStamp(comm, "Init", &t0);
  PetscOptionsBegin(comm, "", "h2mole.c options", "none");
  PetscOptionsGetString(NULL, "-in_dir", in_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-out_dir", out_dir, 100, NULL);
  PetscOptionsGetString(NULL, "-guess_type", guess_type, 10, NULL);
  PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  PetscOptionsGetInt(NULL, "-qmax", &qmax, NULL);
  PetscOptionsGetString(NULL, "-eri", eri_option, 10, NULL);
  ierr = FEMInfCreateFromOptions(&fem, comm); CHKERRQ(ierr);
  PetscOptionsEnd();

  PetscBool s_is_id; FEMInfGetOverlapIsId(fem, &s_is_id);

  Mat s_r1, s_y2, S;
  PrintTimeStamp(comm, "SMat", NULL);
  ierr = FEMInfSetSR1Mat(fem, &s_r1); CHKERRQ(ierr);
  ierr = MatSetDirFile(in_dir, "s_y2mat.dat", &s_y2);   CHKERRQ(ierr);
  if(s_is_id) 
    S = NULL;
  else {
    ierr = MatSetSynthesize3Fast(s_r1, s_r1, s_y2, comm, &S);
  }

  // D2 
  Mat d2_r1, H, D;
  PrintTimeStamp(PETSC_COMM_SELF, "D2Mat", NULL);
  ierr = FEMInfSetD2R1Mat(fem, &d2_r1); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = MatScale(d2_r1, -0.5); CHKERRQ(ierr);
  ierr = MatSetSynthesize3Fast(s_r1, d2_r1, s_y2, comm, &H); CHKERRQ(ierr);
  ierr = MatSetSynthesize3Fast(d2_r1, s_r1, s_y2, comm, &D); CHKERRQ(ierr);
  MatAXPY(H, 1.0, D, DIFFERENT_NONZERO_PATTERN);
  MatDestroy(&d2_r1); MatDestroy(&D);

  // L
  Mat l_r1; 
  PrintTimeStamp(PETSC_COMM_SELF, "LMat", NULL);
  ierr = FEMInfSetR2invR1Mat(fem, &l_r1); CHKERRQ(ierr);

  char l1_path[100]; sprintf(l1_path, "%s/l_1_y2mat.dat", in_dir); 
  FILE *fl1 = fopen(l1_path, "r");
  if(fl1 != NULL) {
    Mat l_1_y2, L;
    ierr = MatCreateFromCOOFormatFileHandler(fl1, &l_1_y2); CHKERRQ(ierr);
    ierr = MatSetSynthesize3Fast(l_r1, s_r1, l_1_y2, comm, &L);
    MatAXPY(H, 0.5, L, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&L); MatDestroy(&l_1_y2); fclose(fl1);
  }

  char l2_path[100]; sprintf(l2_path, "%s/l_2_y2mat.dat", in_dir); 
  FILE *fl2 = fopen(l2_path, "r");
  if(fl2 != NULL) {
    Mat l_2_y2, L;
    ierr = MatCreateFromCOOFormatFileHandler(fl1, &l_2_y2); CHKERRQ(ierr);
    ierr = MatSetSynthesize3Fast(s_r1, l_r1, l_2_y2, comm, &L);
    MatAXPY(H, 0.5, L, DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&L); MatDestroy(&l_2_y2); fclose(fl2);
  }
  MatDestroy(&l_r1);

  // e-n
  PrintTimeStamp(PETSC_COMM_SELF, "ENMat", NULL);
  PetscScalar a = bond_length/2.0;
  for(int q = 0; q < qmax; q++) {
    char path1[100]; char path2[100]; 
    sprintf(path1, "%s/p%d_A1_y2mat.dat", in_dir, q);
    sprintf(path2, "%s/p%d_A2_y2mat.dat", in_dir, q);
    FILE* f1 = fopen(path1, "r"); FILE* f2 = fopen(path2, "r");
    if(f1 != NULL && f2 != NULL) {
      PetscPrintf(comm, "%d ", q);
      Mat q_r1;
      ierr = FEMInfSetENR1Mat(fem, q, a, &q_r1);  CHKERRQ(ierr);

      Mat pq_1_y2;
      ierr = MatCreateFromCOOFormatFileHandler(f1, &pq_1_y2); CHKERRQ(ierr);

      Mat V1;
      ierr = MatSetSynthesize3Fast(q_r1, s_r1, pq_1_y2, comm, &V1); CHKERRQ(ierr);
      MatAXPY(H, -2.0, V1, DIFFERENT_NONZERO_PATTERN);
      MatDestroy(&pq_1_y2); MatDestroy(&V1);

      Mat pq_2_y2;
      ierr = MatCreateFromCOOFormatFileHandler(f2, &pq_2_y2); CHKERRQ(ierr);

      Mat V2;
      ierr = MatSetSynthesize3Fast(s_r1, q_r1, pq_2_y2, comm, &V2); CHKERRQ(ierr);
      MatAXPY(H, -2.0, V2, DIFFERENT_NONZERO_PATTERN);
      MatDestroy(&pq_2_y2); MatDestroy(&V2); 

      MatDestroy(&q_r1); fclose(f1); fclose(f2);
    }
  }
  PetscPrintf(comm, "\n");

  // e-e
  if(strcmp(eri_option, "direct") == 0) {
    PrintTimeStamp(comm, "EEMat", NULL);    
    for(int q = 0; q < qmax; q++) {
      char path[100]; sprintf(path, "%s/p%d_12_y2mat.dat", in_dir, q);
      FILE *f = fopen(path, "r"); 
      if(f != NULL) {
	PetscPrintf(comm, "%d ", q);
	Mat r2, y2, V; 
	ierr = FEMInfSetEER2Mat(fem, q, &r2); CHKERRQ(ierr);
	ierr = MatCreateFromCOOFormatFileHandler(f, &y2); CHKERRQ(ierr);
	ierr = MatSetSynthesizeFast(r2, y2, comm, &V); CHKERRQ(ierr);
	ierr = MatAXPY(H, 1.0, V, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	MatDestroy(&r2); MatDestroy(&y2); MatDestroy(&V);
	fclose(f);
      }
    }
    PetscPrintf(comm, "\n");
  }

  // Solve
  EPS eps; 
  time_t t_solve; PrintTimeStamp(PETSC_COMM_SELF, "eps", &t_solve);
  EPSCreate(comm, &eps); EPSSetOperators(eps, H, S);
  EPSSetTarget(eps, -4.0); 
  if(s_is_id) 
    EPSSetProblemType(eps, EPS_HEP); 
  else
    EPSSetProblemType(eps, EPS_GHEP); 

  if(strcmp(guess_type, "vec") == 0) {
    PrintTimeStamp(comm, "guess", NULL);
    PetscViewer viewer;
    Vec *guess; guess = (Vec*)malloc(sizeof(Vec)*1);
    VecCreate(comm, &guess[0]);
    char path[100]; sprintf(path, "%s/guess.vec.dat", in_dir);
    PetscViewerBinaryOpen(comm, path, FILE_MODE_READ, &viewer);
    ierr = VecLoad(guess[0], viewer); CHKERRQ(ierr);    
    EPSSetInitialSpace(eps, 1, guess);
    VecDestroy(&guess[0]);
    free(guess);
  }
  EPSSetFromOptions(eps);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE); 
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  // Output
  time_t t1; PrintTimeStamp(PETSC_COMM_SELF, "Output", &t1);
  FEMInfFPrintf(fem, stdout, 0);
  PetscPrintf(comm, "\n==== Condition ====\n");
  PetscPrintf(comm, "in_dir: %s\n", in_dir);
  PetscPrintf(comm, "out_dir: %s\n", out_dir);
  PetscPrintf(comm, "bond_length: %f\n", bond_length);
  PetscPrintf(comm, "qmax: %d\n", qmax);
  PetscPrintf(comm, "ERI: %s\n", eri_option);
  
  PetscPrintf(comm, "\n==== Time ====\n");
  PetscPrintf(comm, "t(mat): %f\n", (double)(t_solve-t0));
  PetscPrintf(comm, "t(diag): %f\n", (double)(t1-t_solve));

  EPSType eps_type; EPSGetType(eps, &eps_type);
  PetscPrintf(comm, "eps_type: %s\n", eps_type);

  PetscPrintf(comm, "\n==== Results ====\n");
  int nconv;
  PetscScalar kr;
  Vec xr;
  ierr = MatCreateVecs(H, NULL, &xr); CHKERRQ(ierr);
  EPSGetConverged(eps, &nconv);
  for(int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &kr, NULL, xr, NULL);
    PetscPrintf(comm, "eig%i: %f\n", i, kr);
    PetscViewer viewer;
    char path[100]; sprintf(path, "%s/eig%d.vec.dat", out_dir, i);
    PetscViewerBinaryOpen(comm, path, FILE_MODE_WRITE, &viewer);
    VecView(xr, viewer);
  }  

  // Destroy
  VecDestroy(&xr);
  EPSDestroy(&eps);
  MatDestroy(&H); 
  if(S!=NULL)
    MatDestroy(&S);
  FEMInfDestroy(&fem);

  return 0;
}
*/
