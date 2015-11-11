#include <rescol/oce1.h>
#include <rescol/viewerfunc.h>
#include <rescol/eeps.h>

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

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = PETSC_COMM_SELF;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  PrintTimeStamp(comm, "Init", NULL);
  OCE1 oce1; OCE1Create(comm, &oce1); 
  EEPS eps;  EEPSCreate(comm, &eps);
  double bond_length = 2.0;
  PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
  PetscViewerFormat format;
  ViewerFunc viewer_func; ViewerFuncCreate(comm, &viewer_func);
  
  PetscOptionsBegin(PETSC_COMM_SELF, "", "h2plus.c options", "none");
  ierr = PetscOptionsGetReal(NULL, "-bond_length", &bond_length, NULL);
  ierr = OCE1SetFromOptions(oce1); CHKERRQ(ierr);
  ierr = EEPSSetFromOptions(eps); CHKERRQ(ierr);
  ierr = PetscOptionsGetViewer(comm, NULL, "-view", &viewer, &format, NULL);
  CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(viewer_func); CHKERRQ(ierr);
  
  PetscOptionsEnd();  

  // other conf
  PetscBool s_is_id;

  // Matrix
  PrintTimeStamp(comm, "Mat", NULL);
  Mat H; OCE1CreateMat(oce1, &H);
  Mat S; OCE1CreateMat(oce1, &S); 
  ierr = OCE1TMat(oce1, H); CHKERRQ(ierr);
  ierr = OCE1PlusVneMat(oce1, bond_length/2.0, 1.0, H); CHKERRQ(ierr);
  OCE1SMat(oce1, S, &s_is_id);

  // Solve
  PrintTimeStamp(comm, "EPS", NULL);
  if(s_is_id)
    EEPSSetOperators(eps, H, NULL);
  else
    EEPSSetOperators(eps, H, S);
  EEPSSolve(eps);  

  // Write
  if(ViewerFuncIsActive(viewer_func)) {
    Vec c; MatCreateVecs(H, &c, NULL);
    ierr = EPSGetEigenpair(eps->eps, 0, NULL, NULL, c, NULL); CHKERRQ(ierr);
    ierr = OCE1ViewFunc(oce1, c, viewer_func); CHKERRQ(ierr);
    ierr = ViewerFuncView(viewer_func, viewer);
  }

  // output
  PrintTimeStamp(comm, "Output", NULL);
  OCE1View(oce1, viewer);
  return 0;
}

