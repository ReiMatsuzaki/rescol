
#include <rescol/angmoment.h>


// --- Utils ---
PetscReal wigner3j(int a, int b, int c, int d, int e, int f) {
  return gsl_sf_coupling_3j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}
PetscReal wigner6j(int a, int b, int c, int d, int e, int f) {
  return gsl_sf_coupling_6j(2*a, 2*b, 2*c, 2*d, 2*e, 2*f);
}
PetscBool TriangleQ(int j1, int j2, int j3) {
  return (
	  (abs(j1-j2) <= j3 && j3 <= j1+j2) ||
	  (abs(j2-j3) <= j1 && j1 <= j2+j3) ||
	  (abs(j3-j1) <= j2 && j2 <= j3+j1)
	  );
  return 0;
}
int sign(int n) {
  if(n%2==0)
    return 1;
  else
    return -1;
}
