#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include <vector>
#include <boost/numpy.hpp>

/*
  Bachau, H., Cormier, E., Decleva, J., Hansen, J. E. and Martin, F.
  "Applications of B -splines in atomic and molecular"
  Reports on Progree in Physics *64* (2001) 1815

  Warning:
  almost all variable are same as the above paper,
  but index starts 0 in this program.
 */
double CalcBSpline(int order, double* ts, int i, double x);
double CalcDerivBSpline(int order, double* ts, int i, double x);
double ERI_ele(double* vs, int num_basis, int num_quad, 
	       int i0, int i1, int j0, int j1,
	       int a, int b, int c, int d, 
	       double* ws, double* rij);

#endif
