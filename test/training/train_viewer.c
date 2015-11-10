#include <slepceps.h> 

static char help[] = "solve one particle eigen energy problem";

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  PetscViewer viewer;
  PetscViewerFormat format;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  PetscOptionsGetViewer(comm, NULL, "-my_view", &viewer, &format, NULL);

  if(format == PETSC_VIEWER_ASCII_INFO)
    PetscPrintf(comm, "ascii_info is fountd\n");
  else if(format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    PetscPrintf(comm, "ascii_info_detail is fountd\n");
  } else
    PetscPrintf(comm, "other");

  PetscViewerType type;
  PetscViewerGetType(viewer, &type);
  PetscPrintf(comm, "type: %s\n", type);

  if(strcmp(type, "ascii") == 0) {

    FILE *fp;
    PetscViewerASCIIGetPointer(viewer, &fp);
    PetscFPrintf(comm, fp, "ABC\n");
    PetscFPrintf(comm, fp, "ABC\n");

  } else {
    PetscPrintf(comm, "only ascii is supported\n");
  }
  
  PetscViewerDestroy(&viewer);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
