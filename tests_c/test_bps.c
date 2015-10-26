#include <slepceps.h>
#include "unittest.h"
#include "../src_c/bps.h"

static char help[] = "Unit test for dvr.c \n\n";

PetscErrorCode testLine() {

  PetscErrorCode ierr;
  BPS bps;
  ierr = BPSCreate(&bps, PETSC_COMM_SELF); CHKERRQ(ierr);
  ierr = BPSSetLine(bps, 5.0, 6); CHKERRQ(ierr);

  double *zs;
  int num;
  ierr = BPSGetZs(bps, &zs, &num);

  ASSERT_EQ(6, num)
  
  for(int i = 0; i < num; i++)
    ASSERT_DOUBLE_EQ(i*1.0, zs[i]);

  BPSDestroy(&bps);

  return 0;

}

PetscErrorCode testExp() {

  PetscErrorCode ierr;
  BPS bps;
  ierr = BPSCreate(&bps, PETSC_COMM_SELF); CHKERRQ(ierr);
  ierr = BPSSetFromOptions(bps); CHKERRQ(ierr);

  double *zs;
  int num;
  ierr = BPSGetZs(bps, &zs, &num);

  ASSERT_EQ(6, num);

  BPSFPrintf(bps, stdout, 0);

  BPSDestroy(&bps);

  return 0;

}

int main(int argc, char **args) {
  
  PetscErrorCode ierr;

  SlepcInitialize(&argc, &args, (char*)0, help);
  ierr = testLine(); CHKERRQ(ierr);
  ierr = testExp(); CHKERRQ(ierr);
  SlepcFinalize();
  return 0;
}

