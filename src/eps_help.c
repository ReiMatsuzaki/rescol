#include <slepceps.h>

static char help[] = "for see eps help message";


int main(int argc, char **args) {
  SlepcInitialize(&argc, &args, (char*)0, help); 

  EPS eps;
  EPSCreate(MPI_COMM_SELF, &eps);
  EPSSetFromOptions(eps);
  EPSDestroy(&eps);
  SlepcFinalize();
  return 0;
}

