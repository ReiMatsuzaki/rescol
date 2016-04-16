#include <rescol/oce1.h>
#include <rescol/h2plus.h>
#include <rescol/eeps.h>
#include "unittest.h"

static char help[] = "Unit test for h2plus.c \n\n";

int TestH2plus() {

  MPI_Comm comm = PETSC_COMM_SELF;
  BPS bps; BPSCreate(comm, &bps); BPSSetLine(bps, 100.0, 101);
  BSS bss; BSSCreate(comm, &bss); 
  BSSSetKnots(bss, 5, bps); BSSSetUp(bss);
  FEMInf fem; FEMInfCreate(comm, &fem); FEMInfSetBSS(fem, bss);
  Y1s y1s;  Y1sCreate(comm, &y1s); Y1sSet(y1s, SIGMA, GERADE, 4);
  OCE1 oce; OCE1Create(comm, &oce); OCE1Set(oce, fem, y1s);      

  OceH2plus h2plus; OceH2plusCreate(comm, &h2plus); 
  Mat H, S; 
  PetscBool is_id;
  OceH2plusCalc(h2plus, &H, &S, &is_id);
  
  EEPS eps; EEPSCreate(comm, &eps); EEPSSetOperators(eps, H, S);
  EEpsSetTarget(eps, -2.0);
  

  OCE1Destroy(&oce);
  CtxOceH2PlusDestroy(&h2plus);

  return 0;
}

int main (int argc, char **args) {
  PetscInitialize(&argc, &args, (char*)0, help);
  TestH2plus();
  PetscFinalize();
  return 0;
}

