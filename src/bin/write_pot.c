#include <rescol/pot.h>
#include <rescol/writer.h>

static char help[] = "write potential curve";

int main(int argc, char **args) {
  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  WFWriter writer;
  POT pot;
  char out[100] = "pot.dat";
  int J = 0;
  PetscReal mu = 1.0;

  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  ierr = WFWriterCreateFromOptions(&writer, comm); CHKERRQ(ierr);
  ierr = POTCreateFromOptions(&pot, comm); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL, "-out", out, 100, NULL);
  ierr = PetscOptionsGetInt(NULL, "-J", &J, NULL);
  ierr = PetscOptionsGetReal(NULL, "-mu", &mu, NULL);

  ierr = WFWriterWriteFilePOT(writer, out, pot, J, mu);

  PetscPrintf(comm, "J=%d", J);
  PetscPrintf(comm, "mu=%f", mu);
  WFWriterView(writer);
  POTView(pot); 

  WFWriterDestroy(&writer);
  POTDestroy(&pot);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}

