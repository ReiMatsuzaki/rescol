#include "../include/oce1.h"
#include "../include/viewerfunc.h"
#include "../include/eeps.h"

static char help[] = "solve H2^+ problem";
/*
  -eps_nev : # of necessary eigen pairs
  -eps_ncv : # of column vectors to be used by the solutions
  -eps_max_it : maximum iteration 
  -eps_mpd : maximum dim of projected problem

  -bond_length : 
  -lmax :
  -gerade or -ungerade

  -fem_type
  -fd_num  -fd_xmax
  -bss_order 
  -dvr_nq
  -bps_type -bps_num_zs -bps_zmax
*/

typedef struct {
  MPI_Comm comm;
  PetscReal R;    // bond length
  OCE1      oce;  // grid mehtod for one center expansion
  ViewerFunc viewer_func;
  EEPS        eps;
  char path_in[100];
  char path_out[100];
} p_EigH2plus;
typedef p_EigH2plus *EigH2plus;

PetscErrorCode EigH2plusCreate(MPI_Comm comm, EigH2plus *p_self) {

  PetscErrorCode ierr;
  PrintTimeStamp(comm, "Init", NULL);

  EigH2plus self;
  ierr = PetscNew(&self); CHKERRQ(ierr);
  self->comm = comm;
  *p_self = self;

  self->R = 2.0;
  ierr = OCE1Create(comm, &self->oce); CHKERRQ(ierr);
  ierr = EEPSCreate(comm, &self->eps); CHKERRQ(ierr);

  return 0;
}
PetscErrorCode EigH2plusSetFromOptions(EigH2plus self) {
  PetscErrorCode ierr;
  PetscBool find = PETSC_TRUE;
  PrintTimeStamp(self->comm, "Set", NULL);
  PetscOptionsBegin(self->comm, "", "eig_h2plus.c options", "none");

  // -- set bond length --
  ierr = PetscOptionsGetReal(NULL, NULL, "-R", &self->R, &find); CHKERRQ(ierr);
  if(!find) {
    SETERRQ(self->comm, 1, "-R not found.");
  }

  // -- set grid method --
  ierr = OCE1SetFromOptions(self->oce); CHKERRQ(ierr);
  // -- set function viewer --
  //  ierr = ViewerFuncSetFromOptions(self->viewer_func, &find); CHKERRQ(ierr);
  //  self->use_func_view = find;

  // -- EPS --
  //  ierr = EEPSSetTarget(self->eps, -1.2);  CHKERRQ(ierr);
  ierr = EEPSSetFromOptions(self->eps); CHKERRQ(ierr);
  
  // -- in/out --
  ierr = PetscOptionsGetString(NULL, NULL, "-in", self->path_in, 100, &find); CHKERRQ(ierr);
  if(!find) {
    SETERRQ(self->comm, 1, "-in not found");
  }
  PetscOptionsGetString(NULL, NULL, "-out", self->path_out,
			100, &find);
  if(!find) {
    SETERRQ(self->comm, 1, "-out not found");
  }  
  PetscOptionsEnd();

  return 0;
}
PetscErrorCode EigH2plusPrintIn(EigH2plus self, PetscViewer v) {
  
  PrintTimeStamp(self->comm, "PrintIn", NULL);

  PetscViewerASCIIPrintf(v, "R:%f\n", self->R);

  PetscViewerASCIIPrintf(v, "oce: ");
  OCE1View(self->oce, v);

  //  PetscViewerASCIIPrintf(v, "viewer_func: ");
  //  ViewerFuncView(self->viewer_func, v);

  PetscPrintf(self->comm, "in: %s\n", self->path_in);
  PetscPrintf(self->comm, "out: %s\n", self->path_out);
  return 0;
}
PetscErrorCode EigH2plusCalc(EigH2plus self, PetscViewer v) {

  PetscErrorCode ierr;  

  PetscBool s_is_id;

  // read initial guess
  PrintTimeStamp(self->comm, "read in", NULL);
  Vec c_guess;
  PetscViewer v_in;
  ierr = PetscViewerBinaryOpen(self->comm, self->path_in, FILE_MODE_READ, &v_in); CHKERRQ(ierr);  
  ierr = OCE1CreateVec(self->oce, &c_guess); CHKERRQ(ierr);
  ierr = VecLoad(c_guess, v_in); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&v_in); CHKERRQ(ierr);

  // Matrix
  PrintTimeStamp(self->comm, "Mat", NULL);
  Mat H;
  Mat S;
  ierr = OCE1TMat(self->oce, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);
  ierr = OCE1PlusVneMat(self->oce, self->R/2.0, 1.0, H); CHKERRQ(ierr);
  ierr = OCE1SMat(self->oce, MAT_INITIAL_MATRIX, &S, &s_is_id); CHKERRQ(ierr);


  // Solve
  PrintTimeStamp(self->comm, "EPS", NULL);
  ierr = EPSSetInitialSpace(self->eps->eps, 1, &c_guess); CHKERRQ(ierr);
  PrintTimeStamp(self->comm, "EPS", NULL);
  if(s_is_id)
    EEPSSetOperators(self->eps, H, NULL);
  else
    EEPSSetOperators(self->eps, H, S);
  ierr = EEPSSolve(self->eps); CHKERRQ(ierr);
  
  // Write
  PrintTimeStamp(self->comm, "Write", NULL);
  /*
  if(ViewerFuncIsActive(self->viewer_func)) {
    Vec c; MatCreateVecs(H, &c, NULL);
    ierr = EPSGetEigenpair(self->eps->eps, 0, NULL, NULL, c, NULL); CHKERRQ(ierr);
    ierr = OCE1ViewFunc(self->oce, c, self->viewer_func); CHKERRQ(ierr);

  }
  */
  return 0;
}
PetscErrorCode EigH2plusPrintOut(EigH2plus self, PetscViewer v) {

  int nconv;
  EPSGetConverged(self->eps->eps, &nconv);
  if(nconv == 0) {
    PetscViewerASCIIPrintf(v, "Failed to compute eigen values \n");
    return 0;
  }

  PetscViewerASCIIPrintf(v, "nconv: %d\n", nconv);
  for(int i = 0; i < nconv; i++) {
    PetscScalar ene;
    Vec c; OCE1CreateVec(self->oce, &c);
    EPSGetEigenpair(self->eps->eps, i, &ene, NULL, c, NULL);
    PetscViewerASCIIPrintf(v, "E%d = (%f, %f)\n", i, creal(ene), cimag(ene));
    VecDestroy(&c);
  }
  
  
  return 0;
}
int main(int argc, char **args) {
  
  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscViewer v = PETSC_VIEWER_STDOUT_SELF;
  PetscViewerFormat format;

  PetscViewerASCIIPrintf(v, "\n");
  PetscViewerASCIIPrintf(v, "eig_h2plus program.\n");
  PetscViewerASCIIPrintf(v, "Solve eigen value problem ob hydrogen ion\n");
  PetscViewerASCIIPrintf(v, "\n");

  EigH2plus h2plus;
  ierr = EigH2plusCreate(comm, &h2plus); CHKERRQ(ierr);
  ierr = EigH2plusSetFromOptions(h2plus); CHKERRQ(ierr);
  ierr = EigH2plusPrintIn(h2plus, v); CHKERRQ(ierr);
  ierr = EigH2plusCalc(h2plus, v); CHKERRQ(ierr);
  ierr = EigH2plusPrintOut(h2plus, v); CHKERRQ(ierr);

  return 0;
}

