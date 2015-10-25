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
