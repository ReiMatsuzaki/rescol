#include <iostream>
#include <stdexcept>
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

np::ndarray ERI(const np::ndarray& xs, int L) {

  // type error check
  if(xs.get_nd() != 1) 
    throw std::runtime_error("xs must be one-dim");
  if(xs.get_dtype() != np::dtype::get_builtin<double>())
    throw std::runtime_error("xs must be double array");
  
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  size_t num_x = xs.shape(0);
  double* dd_ys = new double[num_x*num_x];
  for(int i = 0; i < num_x; i++)
    for(int j = 0; j < num_x; j++) {
      double s, g; 
      if(d_xs[i] < d_xs[j]) {
	s = d_xs[i]; g = d_xs[j];
      } else {
	g = d_xs[i]; s = d_xs[j];
      }
      double acc = 1.0/g;
      double sg = s/g;
      for(int L0 = 0; L0 < L; L0++)
	acc *= sg;
      dd_ys[i+num_x*j] = acc;
    }
  return np::from_data(dd_ys,
		       np::dtype::get_builtin<double>(),
		       bp::make_tuple(num_x, num_x),
		       bp::make_tuple(sizeof(double)),
		       bp::object());
}

double PartialDot(const np::ndarray& xs, const np::ndarray& ys, int i0, int i1) {
  // type error check
  if(xs.get_nd() != 1 || ys.get_nd() != 1) 
    throw std::runtime_error("xs and ys must be one-dim");
  if(xs.get_dtype() != np::dtype::get_builtin<double>() ||
     ys.get_dtype() != np::dtype::get_builtin<double>())
    throw std::runtime_error("xs and ys must be double array");

  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ys = reinterpret_cast<double*>(ys.get_data());

  double res(0.0);
  for(int i = i0; i < i1; i++)
    res += d_xs[i] * d_ys[i];
  return res;
}

double PartialDot3(const np::ndarray& xs, const np::ndarray& ys, 
		   const np::ndarray& zs, int i0, int i1) {
  // type error check
  if(xs.get_nd() != 1 || ys.get_nd() != 1 || zs.get_nd() != 1) 
    throw std::runtime_error("xs and ys must be one-dim");
  if(xs.get_dtype() != np::dtype::get_builtin<double>() ||
     ys.get_dtype() != np::dtype::get_builtin<double>() ||
     zs.get_dtype() != np::dtype::get_builtin<double>())
    throw std::runtime_error("xs and ys must be double array");

  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ys = reinterpret_cast<double*>(ys.get_data());
  double* d_zs = reinterpret_cast<double*>(zs.get_data());

  double res(0.0);
  for(int i = i0; i < i1; i++)
    res += d_xs[i] * d_ys[i] * d_zs[i];
  return res;
}

double Dot_abwAcdw(const np::ndarray& as, const np::ndarray& bs, 
			const np::ndarray& cs, const np::ndarray& ds,
		   const np::ndarray& ws, const np::ndarray& A,
		   int i0, int i1, int j0, int j1) {
  double* d_as = reinterpret_cast<double*>(as.get_data());
  double* d_bs = reinterpret_cast<double*>(bs.get_data());
  double* d_cs = reinterpret_cast<double*>(cs.get_data());
  double* d_ds = reinterpret_cast<double*>(ds.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());  
  double* d_A = reinterpret_cast<double*>(A.get_data());  
  double res(0.0);
  int n = as.shape(0);
  for(int i = i0; i < i1; i++)
    for(int j = j0; j< j1; j++)
      res += d_as[i]*d_cs[i]*d_ws[i]*d_A[i+n*j]*d_bs[j]*d_ds[j]*d_ws[j];
  return res;
}

void Non0QuadIndex(int a, int c, int k, int nq, int* i0, int* i1) {
  *i0 = a<c ? (c-k+2)*k : (a-k+2)*k;
  if(*i0<0)
    *i0 = 0;
  *i1 = a<c ? (a+k-1)*k : (c+k-1)*k;
  if(*i1>nq)
    *i1=nq;
}

bp::tuple Non0QuadIndexPy(int a, int c, int k, int nq) {
  int i0, i1;
  Non0QuadIndex(a, c, k, nq, &i0, &i1);
  return bp::make_tuple(i0, i1);
}

bp::tuple ERI_mat(const np::ndarray& vals, const np::ndarray& xs, 
		  const np::ndarray& ws, int L, int k) {
		  
  double* d_vals = reinterpret_cast<double*>(vals.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());
  int nq = xs.shape(0);
  int nb = vals.shape(0)/nq;
  double* sg_ij = new double[nq*nq];
  for(int i = 0; i < nq; i++)
    for(int j = 0; j < nq; j++) {
      double s, g;
      if (i<j) {
	s = d_xs[i]; g = d_xs[j];
      } else {
	s = d_xs[j]; g = d_xs[i];
      }
      double sg = s/g;
      sg_ij[j+nq*i] = 1.0/g;
      for(int ll = 0; ll < L; ll++)
	sg_ij[j+nq*i] *= sg;
    }

  int num_ele = (nb-k)*(2*k-1)+k*k;
  num_ele = num_ele * num_ele;
  double* data = new double[num_ele];
  int* row = new int[num_ele];
  int* col = new int[num_ele];
  int idx(0);
  for(int a = 0; a < nb; a++) {
    for(int b= 0; b < nb; b++) {
      for(int c = 0; c < nb; c++) {
	for(int d = 0; d < nb; d++) {
	  if(std::abs(a-c) < k && std::abs(b-d) < k) {	    
	    if(idx >= num_ele) 
	      throw std::runtime_error("Exceed index in ERI_mat");
	    int i0, i1, j0, j1;
	    Non0QuadIndex(a, c, k, nq, &i0, &i1);
	    Non0QuadIndex(b, d, k, nq, &j0, &j1);
	    data[idx] = ERI_ele(d_vals, nb, nq, 
				i0, i1, j0, j1,
				a, b, c, d, 
				d_ws, sg_ij);
	    col[idx] = a*nb+b;
	    row[idx] = c*nb+d;
	    idx++;
	  }
	}
      }
    }
  }
  std::cout << "c++, data[0]:" << data[0] << std::endl;
  np::ndarray np_data = np::from_data(data,
			  np::dtype::get_builtin<double>(),
			  bp::make_tuple(num_ele),
			  bp::make_tuple(sizeof(double)),
			  bp::object());
  np::ndarray np_row = np::from_data(row,
			 np::dtype::get_builtin<int>(),
			 bp::make_tuple(num_ele),
			 bp::make_tuple(sizeof(int)),
			 bp::object());
  np::ndarray np_col = np::from_data(col,
			 np::dtype::get_builtin<int>(),
			 bp::make_tuple(num_ele),
			 bp::make_tuple(sizeof(int)),
			 bp::object());
  return bp::make_tuple(np_data, np_row, np_col);
}

BOOST_PYTHON_MODULE(bspline_bind) {
  Py_Initialize();
  np::initialize();
  def( "calc_bspline_xs", CalcBSplineNumpy);
  def("calc_deriv_bspline_xs", CalcDerivBSplineNumpy);
  def("ra_inv", RAInv);
  def("eri", ERI);
  def("pdot", PartialDot);
  def("pdot3", PartialDot3);
  def("dot_abwAcdw", Dot_abwAcdw);
  def("non0_quad_index", Non0QuadIndexPy);
  def("eri_mat", ERI_mat);
}

