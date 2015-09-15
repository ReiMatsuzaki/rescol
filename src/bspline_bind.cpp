#include <iostream>
#include <boost/python.hpp>
#include "bspline.hpp"

namespace {
  using namespace boost::python;
  namespace bp = boost::python;
  namespace np = boost::numpy;
  using namespace std;
}

// bindings for numpy
np::ndarray CalcBSplineNumpy(int order, const np::ndarray& ts,
			     int index, const np::ndarray& xs) {
  double* d_ts = reinterpret_cast<double*>(ts.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  size_t num_x = xs.shape(0);
  double* d_ys = new double[num_x];
  for(int i = 0; i < num_x; i++) 
    d_ys[i] = CalcBSpline(order, d_ts, index, d_xs[i]);
  return np::from_data( d_ys, 
			np::dtype::get_builtin<double>(),
			bp::make_tuple(num_x),
			bp::make_tuple(sizeof(double)),
			bp::object());
}
np::ndarray CalcDerivBSplineNumpy(int order, const np::ndarray& ts,
				  int index, const np::ndarray& xs) {
  double* d_ts = reinterpret_cast<double*>(ts.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  size_t num_x = xs.shape(0);
  double* d_ys = new double[num_x];
  for(int i = 0; i < num_x; i++) 
    d_ys[i] = CalcDerivBSpline(order, d_ts, index, d_xs[i]);
  return np::from_data( d_ys, 
			np::dtype::get_builtin<double>(),
			bp::make_tuple(num_x),
			bp::make_tuple(sizeof(double)),
			bp::object());
}

np::ndarray RAInv(const np::ndarray& xs, int L, double a) {
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  size_t num_x = xs.shape(0);
  double* d_ys = new double[num_x];
  for(int i = 0; i < num_x; i++) {
    double s;
    double g;
    if(a < d_xs[i]) {
      s = a;
      g = d_xs[i];
    } else {
      s = d_xs[i];
      g = a;
    }
    double acc = 1.0/g;
    double sg = s/g;
    for(int L0 = 0; L0 < L; L0++)
      acc *= sg;
    d_ys[i] = acc;
  }
  return np::from_data( d_ys,
			np::dtype::get_builtin<double>(),
			bp::make_tuple(num_x),
			bp::make_tuple(sizeof(double)),
			bp::object());
}

BOOST_PYTHON_MODULE(bspline_bind) {
  Py_Initialize();
  np::initialize();
  def( "calc_bspline_xs", CalcBSplineNumpy);
  def("calc_deriv_bspline_xs", CalcDerivBSplineNumpy);
  def("ra_inv", RAInv);
}

