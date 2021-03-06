#include <cmath>
#include <vector>
#include <stdexcept>
#include "bspline.hpp"

namespace {
  using std::abs;
}

double CalcBSpline(int order, double* ts, int i, double x) {
  
  if(order < 1) {
    std::domain_error err("bspline.cpp, CalcBSpline, order must be positive.");
    throw err;
  }

  if(i < 0) {
    std::domain_error err("bspline.cpp, CalcBSpline, index must be positive");
    throw err;    
  }

  if(x < ts[i] || ts[i+order] < x)
    return 0.0;

  if(order == 1) {
    if(ts[i] <= x && x < ts[i+1])
      return 1.0;
    else
      return 0.0;
  } else {
    double ti    = ts[i];
    double ti1   = ts[i+1];
    double tidm1 = ts[i+order-1];
    double tid   = ts[i+order];
    double bs0 = CalcBSpline(order-1, ts, i,   x);
    double bs1 = CalcBSpline(order-1, ts, i+1, x);
    
    double acc(0.0);
    if(std::abs(bs0) > 0.000000001)
      acc += (x - ti) / (tidm1 - ti) * bs0;
    if(abs(bs1) > 0.000000001)
      acc += (tid - x) / (tid - ti1) * bs1;
    return acc;
  }
}

double CalcDerivBSpline(int order, double* ts, int i, double x) {
  
  if(order < 1) {
    std::domain_error err("bspline.cpp, CalcDerivBSpline, order must be positive.");
    throw err;
  }

  if(i < 0) {
    std::domain_error err("bspline.cpp, CalcDerivBSpline, index must be positive");
    throw err;    
  }
  
  if(order == 1)
    return 0.0;
  else {

    int k(order);
    double bs0 = CalcBSpline(k-1, ts, i, x);
    double bs1 = CalcBSpline(k-1, ts, i+1, x);
    double eps = 0.000000001;
    double acc(0.0);
    if(std::abs(bs0) > eps)
      acc += (k-1)/(ts[i+k-1]-ts[i]) * bs0;
    if(std::abs(bs1) > eps)
      acc -= (k-1)/(ts[i+k]-ts[i+1]) * bs1;
    return acc;
  }
}

double ERI_ele(double* vs, int num_basis, int num_quad,
	       int i0, int i1, int j0, int j1,
	       int a, int b, int c, int d,
	       double* ws, double* rij) {

  double res(0.0);
  for(int i = i0; i < i1; i++)
    for(int j = j0; j < j1; j++) {
      res += ws[i]*ws[j]*rij[i*num_quad+j]*
	vs[a*num_quad+i]*vs[c*num_quad+i]*
	vs[b*num_quad+j]*vs[d*num_quad+j];
    }
  return res;
}
	       


