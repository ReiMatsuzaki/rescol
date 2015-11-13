#include <rescol/fem_inf.h>
#include <rescol/synthesize.h>

static char help[] = "Time comparing for MatMatSynthesize";

int main (int argc, char **args) {
  PetscInitialize(&argc, &args, (char*)0, help);
  MPI_Comm comm = PETSC_COMM_SELF;
  FEMInf fem; FEMInfCreate(comm, &fem);
  FEMInfSetFromOptions(fem);

  clock_t t0 = clock();
  Mat S; FEMInfCreateMat(fem, 1, &S); FEMInfSR1Mat(fem, S);
  Mat V; FEMInfCreateMat(fem, 1, &V); FEMInfENR1Mat(fem, 0, 0.0, V);
  Mat U; MatMatSynthesize(S, V, 1.0, MAT_INITIAL_MATRIX, &U);
  clock_t t1 = clock();
    //  time_t t1; PrintTimeStamp(comm, "End", &t1);
  PetscPrintf(comm, "t=%f\n",(double)(t1-t0)/CLOCKS_PER_SEC);
  PetscFinalize();
  return 0;
}
