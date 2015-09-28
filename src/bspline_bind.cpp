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

/* TO BE REMOVED
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
*/

/* to be removed
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
*/

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

void CalcERI_L(double* xs, int num, int L, double** eri_L) {
  *eri_L = new double[num*num];
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      double s, g;
      if (i<j) {
	s = xs[i]; g = xs[j];
      } else {
	s = xs[j]; g = xs[i];
      }
      double sg = s/g;
      (*eri_L)[j+num*i] = 1.0/g;
      for(int ll = 0; ll < L; ll++)
	(*eri_L)[j+num*i] *= sg;
    }
}

bp::tuple ERI_mat(const np::ndarray& vals, const np::ndarray& xs, 
		  const np::ndarray& ws, int L, int k) {
  double* d_vals = reinterpret_cast<double*>(vals.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());
  int nq = xs.shape(0);
  int nb = vals.shape(0)/nq;
  double* sg_ij;
  CalcERI_L(d_xs, nq, L, &sg_ij);

  int num_ele = (nb-k)*(2*k-1)+k*k;
  num_ele = num_ele * num_ele;

  double* data = new double[num_ele];
  int* row = new int[num_ele];
  int* col = new int[num_ele];
  for (int i = 0; i < num_ele; i++) {
    row[i] = -1; col[i] = -1; 
  }
  double* wvv = new double[nb*nb*nq];
  int* i0s = new int[nb*nb];
  int* i1s = new int[nb*nb];
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int c = c0; c < c1; c++) {
      int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
      i0s[a*nb+c] = i0; i1s[a*nb+c] = i1;
      for(int i = i0; i < i1; i++)
	wvv[a*nb*nq+c*nq+i] = d_ws[i]*d_vals[a*nq+i]*d_vals[c*nq+i];
    }
  }

  int idx(0);
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int c = c0; c < c1; c++) {
      for(int b= 0; b < nb; b++) {
	int d0 = b-k+1; d0 = d0<0? 0:d0;
	int d1 = b+k;   d1 = d1>nb? nb:d1;
	for(int d = d0; d < d1; d++) {

	  if(a >= c && b >= d) {
	    double res(0.0);
	    for(int i = i0s[a*nb+c]; i < i1s[a*nb+c]; i++) {
	      double aci = wvv[a*nb*nq+c*nq+i];
	      for(int j = i0s[b*nb+d]; j < i1s[b*nb+d]; j++) 
		res +=  aci * sg_ij[i*nq+j] * wvv[b*nb*nq+d*nq+j];
	    }
	    
	    if(a > c && b > d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = a*nb+d; row[idx] = c*nb+b; idx++;	  
	      data[idx] = res; col[idx] = c*nb+b; row[idx] = a*nb+d; idx++;
	      data[idx] = res; col[idx] = c*nb+d; row[idx] = a*nb+b; idx++;

	    } else if(a > c && b == d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = c*nb+b; row[idx] = a*nb+d; idx++;

	    } else if (a == c && b > d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = a*nb+d; row[idx] = c*nb+b; idx++;

	    } else if (a == c && b == d ) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	    }
	  }
	}
      }
    }
  }
  if(idx != num_ele) {
    std::cout << "idx does not match num_ele" << std::endl;
    std::cout << "idx: " << idx << std::endl;
    std::cout << "num_ele: " << num_ele << std::endl;
    throw std::runtime_error("idx does not match num_ele");
  }
  delete[] wvv; delete[] i0s;delete[] i1s;
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

bp::tuple ERI_mat8sym(const np::ndarray& vals, const np::ndarray& xs, 
		  const np::ndarray& ws, int L, int k) {
  std::cout << "Eight Symmery" << std::endl;
  double* d_vals = reinterpret_cast<double*>(vals.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());
  int nq = xs.shape(0);
  int nb = vals.shape(0)/nq;
  double* sg_ij;
  CalcERI_L(d_xs, nq, L, &sg_ij);

  int num_ele = (nb-k)*(2*k-1)+k*k;
  num_ele = num_ele * num_ele;

  double* data = new double[num_ele];
  int* row = new int[num_ele];
  int* col = new int[num_ele];
  for (int i = 0; i < num_ele; i++) {
    row[i] = -1; col[i] = -1; 
  }
  double* wvv = new double[nb*nb*nq];
  int* i0s = new int[nb*nb];
  int* i1s = new int[nb*nb];
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int c = c0; c < c1; c++) {
      int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
      i0s[a*nb+c] = i0; i1s[a*nb+c] = i1;
      for(int i = i0; i < i1; i++)
	wvv[a*nb*nq+c*nq+i] = d_ws[i]*d_vals[a*nq+i]*d_vals[c*nq+i];
    }
  }

  int idx(0);
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int c = c0; c < c1; c++) {
      for(int b= 0; b < nb; b++) {
	int d0 = b-k+1; d0 = d0<0? 0:d0;
	int d1 = b+k;   d1 = d1>nb? nb:d1;
	for(int d = d0; d < d1; d++) {

	  int ac = nb*a+c; 
	  int bd = nb*b+d;
	  if(a >= c && b >= d && ac >= bd) {
	    double res(0.0);
	    for(int i = i0s[a*nb+c]; i < i1s[a*nb+c]; i++) {
	      double aci = wvv[a*nb*nq+c*nq+i];
	      for(int j = i0s[b*nb+d]; j < i1s[b*nb+d]; j++) 
		res +=  aci * sg_ij[i*nq+j] * wvv[b*nb*nq+d*nq+j];
	    }
	    
	    if(a > c && b > d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = a*nb+d; row[idx] = c*nb+b; idx++;	  
	      data[idx] = res; col[idx] = c*nb+b; row[idx] = a*nb+d; idx++;
	      data[idx] = res; col[idx] = c*nb+d; row[idx] = a*nb+b; idx++;
	      
	      if(ac > bd) {
		data[idx] = res; col[idx] = b*nb+a; row[idx] = d*nb+c; idx++;
		data[idx] = res; col[idx] = b*nb+c; row[idx] = d*nb+a; idx++;	  
		data[idx] = res; col[idx] = d*nb+a; row[idx] = b*nb+c; idx++;
		data[idx] = res; col[idx] = d*nb+c; row[idx] = b*nb+a; idx++;
	      }
	    } else if(a > c && b == d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = c*nb+b; row[idx] = a*nb+d; idx++;

	      if(ac > bd) {
		data[idx] = res; col[idx] = b*nb+a; row[idx] = d*nb+c; idx++;
		data[idx] = res; col[idx] = d*nb+a; row[idx] = b*nb+c; idx++;
	      }

	    } else if (a == c && b > d) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      data[idx] = res; col[idx] = a*nb+d; row[idx] = c*nb+b; idx++;

	      if(ac > bd) {
		data[idx] = res; col[idx] = b*nb+a; row[idx] = d*nb+c; idx++;
		data[idx] = res; col[idx] = b*nb+c; row[idx] = d*nb+a; idx++;
	      }
	    } else if (a == c && b == d ) {
	      data[idx] = res; col[idx] = a*nb+b; row[idx] = c*nb+d; idx++;
	      if(ac > bd) 
		data[idx] = res; col[idx] = c*nb+d; row[idx] = a*nb+c; idx++;
	    }
	  }
	}
      }
    }
  }
  if(idx != num_ele) {
    std::cout << "idx does not match num_ele" << std::endl;
    std::cout << "idx: " << idx << std::endl;
    std::cout << "num_ele: " << num_ele << std::endl;
    throw std::runtime_error("idx does not match num_ele");
  }
  delete[] wvv; delete[] i0s;delete[] i1s;
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

/* TO BE REMOVED
bp::tuple ERI_matold(const np::ndarray& vals, const np::ndarray& xs, 
		  const np::ndarray& ws, int L, int k) {
  std::cout << "old" << std::endl;  
  double* d_vals = reinterpret_cast<double*>(vals.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());
  int nq = xs.shape(0);
  int nb = vals.shape(0)/nq;
  double* sg_ij;
  CalcERI_L(d_xs, nq, L, &sg_ij);

  int num_ele = (nb-k)*(2*k-1)+k*k;
  num_ele = num_ele * num_ele;

  double* data = new double[num_ele];
  int* row = new int[num_ele];
  int* col = new int[num_ele];
  int idx(0);
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int b= 0; b < nb; b++) {
      int d0 = b-k+1; d0 = d0<0? 0:d0;
      int d1 = b+k;   d1 = d1>nb? nb:d1;
      for(int c = c0; c < c1; c++) {
	int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
	for(int d = d0; d < d1; d++) {
	  int j0, j1; Non0QuadIndex(b, d, k, nq, &j0, &j1);

	  double res(0.0);
	  for(int i = i0; i < i1; i++)
	    for(int j = j0; j < j1; j++) {
	      res += d_ws[i]*d_ws[j]*sg_ij[i*nq+j]*
		d_vals[a*nq+i]*d_vals[c*nq+i]*
		d_vals[b*nq+j]*d_vals[d*nq+j];
	    }
	  data[idx] = res;
	  col[idx] = a*nb+b;
	  row[idx] = c*nb+d;
	  idx++;
	}
      }
    }
  }
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


bp::tuple ERI_matOld(const np::ndarray& vals, const np::ndarray& xs, 
		     const np::ndarray& ws, int L, int k) {
		  
  double* d_vals = reinterpret_cast<double*>(vals.get_data());
  double* d_xs = reinterpret_cast<double*>(xs.get_data());
  double* d_ws = reinterpret_cast<double*>(ws.get_data());
  int nq = xs.shape(0);
  int nb = vals.shape(0)/nq;
  double* sg_ij;
  CalcERI_L(d_xs, nq, L, &sg_ij);

  int num_ele = (nb-k)*(2*k-1)+k*k;
  num_ele = num_ele * num_ele;

  double* data = new double[num_ele];
  int* row = new int[num_ele];
  int* col = new int[num_ele];
  int idx(0);
  for(int a = 0; a < nb; a++) {
    int c0 = a-k+1;   c0 = c0<0? 0 : c0;
    int c1 = a+k;     c1 = c1>nb? nb : c1;
    for(int b= 0; b < nb; b++) {
      int d0 = b-k+1; d0 = d0<0? 0:d0;
      int d1 = b+k;   d1 = d1>nb? nb:d1;
      for(int c = c0; c < c1; c++) {
	int i0, i1; Non0QuadIndex(a, c, k, nq, &i0, &i1);
	for(int d = d0; d < d1; d++) {
	  int j0, j1; Non0QuadIndex(b, d, k, nq, &j0, &j1);

	  double res(0.0);
	  for(int i = i0; i < i1; i++)
	    for(int j = j0; j < j1; j++) {
	      res += d_ws[i]*d_ws[j]*sg_ij[i*nq+j]*
		d_vals[a*nq+i]*d_vals[c*nq+i]*
		d_vals[b*nq+j]*d_vals[d*nq+j];
	    }
	  data[idx] = res;
	  col[idx] = a*nb+b;
	  row[idx] = c*nb+d;
	  idx++;
	}
      }
    }
  }
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
*/

template<class T>
np::ndarray OuterSumTemplate(np::ndarray& a, np::ndarray& b) {

  // type check
  if(np::dtype::get_builtin<T>() != a.get_dtype() ||
     np::dtype::get_builtin<T>() != b.get_dtype()) 
    throw runtime_error("invalid data type for a or b");

  // dimension check
  if(a.get_nd() != 1 || b.get_nd() != 1)
    throw runtime_error("a and b must be 1-dim array");

  T* da = reinterpret_cast<T*>(a.get_data());
  T* db = reinterpret_cast<T*>(b.get_data());
  int na = a.shape(0);
  int nb = b.shape(0);

  int ny = na * nb;
  T* y = new T[ny];
  for(int i = 0; i < na; i++)
    for(int j = 0; j < nb; j++)
      y[nb*i+j] = (int)(da[i]+db[j]);
  np::ndarray res = np::from_data(y,
				  np::dtype::get_builtin<T>(),
				  bp::make_tuple(ny),
				  bp::make_tuple(sizeof(T)),
				  bp::object());
  return res;
}

np::ndarray OuterSum(np::ndarray& a, np::ndarray& b){
  if(np::dtype::get_builtin<int>() == a.get_dtype() &&
     np::dtype::get_builtin<int>() == b.get_dtype()) 
    return OuterSumTemplate<int>(a, b);
  else if(np::dtype::get_builtin<long long>() == a.get_dtype() &&
     np::dtype::get_builtin<long long>() == b.get_dtype()) 
    return OuterSumTemplate<long long>(a, b);
  else 
    throw std::runtime_error("data type for a or b are not supported");
}

/* TO BE REMOVED
np::ndarray OuterSumInt(np::ndarray& a, np::ndarray& b) {
  long long* da = reinterpret_cast<long long*>(a.get_data());
  int na = a.shape(0);
  long long* db = reinterpret_cast<long long*>(b.get_data());
  int nb = b.shape(0);

  int ny = na*nb;
  int* ys = new int[ny];
  for (int i = 0; i < na; i++) 
    for (int j = 0; j < nb; j++){
      ys[nb*i+j] = (int)(da[i]+db[j]); 
  }

  bp::object own;
  np::ndarray res =  np::from_data(ys,
				   np::dtype::get_builtin<int>(),
				   bp::make_tuple(ny),
				   bp::make_tuple(sizeof(int)),
				   own);

  return res;
}
*/

BOOST_PYTHON_MODULE(bspline_bind) {
  Py_Initialize();
  np::initialize();
  
  def( "calc_bspline_xs", CalcBSplineNumpy);
  def("calc_deriv_bspline_xs", CalcDerivBSplineNumpy);
  def("ra_inv", RAInv);
  //  def("eri", ERI);
  //  def("pdot", PartialDot);
  //  def("pdot3", PartialDot3);
  def("dot_abwAcdw", Dot_abwAcdw);
  def("non0_quad_index", Non0QuadIndexPy);
  def("eri_mat", ERI_mat);
  def("outer_sum", OuterSum);
}

