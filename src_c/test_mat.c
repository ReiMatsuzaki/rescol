#include "mat.h"

static char help[] = "Test for mat.c";

int main (int argc, char **args)
{

  PetscErrorCode ierr;
  char h_path[256];
  Mat H;

  // Initialization
  ierr = SlepcInitialize(&argc, &args, (char*)0, help);
  
  //  MatCreate(PETSC_COMM_WORLD, &H);
  PetscOptionsString("-h-mat", "path for H matrix data file", "", h_path, h_path, 256, NULL);

  ierr = MatCreateFromCOOFormatFile(h_path, &H); CHKERRQ(ierr);
  
  ierr = MatView(H, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

  ierr = MatDestroy(&H); CHKERRQ(ierr);

  ierr = SlepcFinalize();
  return 0;

}
