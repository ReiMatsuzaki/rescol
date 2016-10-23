#include "../include/oce1.h"
#include "../include/viewerfunc.h"

static char help[] = "Fit L=0 function in OCE1";

int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);
  
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscBool find;
  PetscViewer v = PETSC_VIEWER_STDOUT_SELF;

  
  PetscPrintf(comm, "\n");
  PetscPrintf(comm, ">>>> fit_oce1 program >>>>\n");
  PetscPrintf(comm, "Fit L=0 radial function in OCE1\n");
  
  OCE1 oce;
  Pot pot;
  KSP ksp;
  Vec c;
  char path_out[100];
  PetscViewerFormat format;
  PetscBool set;

  // -- create --
  PrintTimeStamp(comm, "Init", NULL);
  ierr = OCE1Create(comm, &oce); CHKERRQ(ierr);
  ierr = PotCreate(comm,  &pot); CHKERRQ(ierr);
  ierr = KSPCreate(comm, &ksp);  CHKERRQ(ierr);
  
  // -- read options --
  PrintTimeStamp(comm, "Set", NULL);
  PetscOptionsBegin(comm, "", "fit_oce1.c options", "none");
  ierr = OCE1SetFromOptions(oce); CHKERRQ(ierr);
  ierr = PotSetFromOptions2(pot, "v_", &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, NULL, "-out", path_out,
			       100, &set); CHKERRQ(ierr);
  CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  // -- input error --
  if(pot == NULL) {
    SETERRQ(comm, 1, "-v_pot option is necessary");
  }
  
  // -- print in --
  PrintTimeStamp(comm, "PrintIn", NULL);
  ierr = PetscPrintf(comm, "OCE1: "); CHKERRQ(ierr);  
  ierr = OCE1View(oce, v); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "POT: "); CHKERRQ(ierr);
  ierr = PFView(pot, v); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "out: %s\n", path_out);
  
  // -- calculation --
  PrintTimeStamp(comm, "Calc", NULL);
  ierr = OCE1CreateVec(oce, &c); CHKERRQ(ierr);
  ierr = OCE1Fit(oce, pot, 0, ksp, c); CHKERRQ(ierr);

  // -- write --
  PetscViewer v_out;
  ierr = PetscViewerBinaryOpen(comm, path_out, FILE_MODE_WRITE, &v_out); CHKERRQ(ierr);
  ierr = VecView(c, v_out); CHKERRQ(ierr);
  PetscViewerDestroy(&v_out);

  PetscPrintf(comm, "<<<< fit_oce1 program <<<<\n\n");


  // -- finalize --
  
  return 0;
  
}  
				   
