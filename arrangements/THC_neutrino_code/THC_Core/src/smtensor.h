/*! \file smtensor.h
 *  \author Wolfgang Kastaun
 */

#ifndef SMTENSOR_H
#define SMTENSOR_H

#include <cmath>

namespace whizza {

enum zero_literal {ZERO=0}; //we want to write matrix=ZERO but not matrix= 14.0
enum one_literal {ONE=1};

template<class T, int N> class sm_matrix_sym;
template<class T, int N, bool UP1, bool UP2> class sm_tensor2_sym;

//---------------------------------------------------------------------------------------
//  Raw vector class
//---------------------------------------------------------------------------------------

template<class T, int N> class sm_vector {
  typedef sm_vector<T,N> me_t;
  T v[N];
  public:
  typedef T elem_t;
  enum {SIZE=N};

  sm_vector(){}
  sm_vector(zero_literal z) {zero();}

  T& operator()(int j) {return v[j];}
  const T& operator()(int j) const {return v[j];}

  void assign_prod(const me_t& a,const T& z) {
    for(int i=0;i<N;i++) v[i]=a.v[i]*z;
  }
  void assign_div(const me_t& a, const T& z) {
    assign_prod(a,1.0/z);
  }
  void assign_sum(const me_t& a,const me_t& b) {
    for(int i=0;i<N;i++) v[i]=a.v[i]+b.v[i];
  }
  void assign_diff(const me_t& a,const me_t& b) {
    for(int i=0;i<N;i++) v[i]=a.v[i]-b.v[i];
  }
  void assign_minus(const me_t& a) {
    for(int i=0;i<N;i++) v[i]=-a.v[i];
  }
  T dot(const me_t& a) const {
    T erg=v[0]*a.v[0];
    for (int i=1; i<N; i++) erg += v[i]*a.v[i];
    return erg;
  }
  void assign_prod(const sm_matrix_sym<T, N>& m, const me_t& v);
  void assign_prod(const me_t& v, const sm_matrix_sym<T, N>& m) {assign_prod(m,v);}

  void operator+=(const me_t& a) {assign_sum(*this,a);}
  void operator-=(const me_t& a) {assign_diff(*this,a);}
  void operator*=(const T& z) {assign_prod(*this,z);}
  void operator/=(const T& z) {assign_div(*this,z);}

  me_t operator+(const me_t& a) const {me_t e; e.assign_sum(*this,a); return e;}
  me_t operator-(const me_t& a) const {me_t e; e.assign_diff(*this,a); return e;}
  me_t operator*(const T& z) const {me_t e; e.assign_prod(*this,z); return e;}
  me_t operator/(const T& z) const {me_t e; e.assign_div(*this,z); return e;}

  void zero() {for(int i=0;i<N;i++) v[i]=0.0;}
};

template<class T, int N>
sm_vector<T, N> operator-(sm_vector<T, N> &v)
{
  sm_vector<T, N> erg;
  erg.assign_minus(v);
  return erg;
}

template<class T, int N>
sm_vector<T, N> operator*(const T& a, const sm_vector<T, N> &v1){
  sm_vector<T, N> erg;
  erg.assign_prod(v1,a);
  return erg;
}


//---------------------------------------------------------------------------------------
//  Raw symmetric matrix class
//---------------------------------------------------------------------------------------

template<class T, int N> class sm_matrix_sym {
  typedef sm_matrix_sym<T, N> me_t;
  enum {FLAT_SIZE=(N*(N+1)/2)};

  static int index(int i,int j) {return  j<=i ? j+(i*(i+1))/2 : i+(j*(j+1))/2;}
  sm_vector<T, FLAT_SIZE> c;

  public:
  typedef T elem_t;

  sm_matrix_sym(){}
  sm_matrix_sym(zero_literal z) {c.zero();}
  sm_matrix_sym(one_literal z) {diag(1.0);}

  T& operator()(int i,int j) {return c(index(i,j));}
  const T& operator()(int i,int j) const {return c(index(i,j));}

  void assign_sum(const me_t& a,const me_t& b) {c.assign_sum(a.c, b.c);}
  void assign_diff(const me_t& a,const me_t& b) {c.assign_diff(a.c, b.c);}
  void assign_minus(const me_t& a) {c.assign_minus(a.c);}
  void assign_prod(const me_t& a, const T& z) {c.assign_prod(a.c, z);}
  void assign_div(const me_t& a, const T& z) {c.assign_div(a.c, z);}

  void operator+=(const me_t& a) {assign_sum(*this, a);}
  void operator-=(const me_t& a) {assign_diff(*this, a);}
  void operator*=(const T& z) {assign_prod(*this, z);}
  void operator/=(const T& z) {assign_div(*this, z);}

  me_t operator+(const me_t &a) const {me_t e; e.assign_sum(*this,a); return e;}
  me_t operator-(const me_t &a) const {me_t e; e.assign_diff(*this,a); return e;}
  me_t operator*(const T& z) const {me_t e; e.assign_prod(*this,z); return e;}
  me_t operator/(const T& z) const {me_t e; e.assign_div(*this,z); return e;}
  void diag(const T& d);

  T bilinear(const sm_vector<T, N>& vl, const sm_vector<T, N>& vr) const;
  T bilinear(const sm_vector<T, N>& v) const;
};

template<class T, int N>
void sm_matrix_sym<T, N>::diag(const T& d) {
  for (int i=0;i<N;i++) {
    (*this)(i,i)=d;
    for (int j=0;j<i;j++) (*this)(i,j)=0.0;
  }
}

template<class T,int N> sm_matrix_sym<T, N> operator-(sm_matrix_sym<T, N> &v)
{
  sm_matrix_sym<T, N> erg;
  erg.assign_minus(v);
  return erg;
}

template<class T, int N>
sm_matrix_sym<T, N> operator*(const T& a, const sm_matrix_sym<T, N>& m){
  sm_matrix_sym<T, N> erg;
  erg.assign_prod(m, a);
  return erg;
}

//---------------------------------------------------------------------------------------
// Vector times matrix
//---------------------------------------------------------------------------------------

template<class T, int N>
void sm_vector<T, N>::assign_prod(const sm_matrix_sym<T, N>& m, const me_t& w)
{
  for (int i=0; i<N; i++) {
    v[i] = w(0)*m(i,0);
    for (int j=1; j<N; j++) v[i] += w(j) * m(i,j);
  }
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_matrix_sym<T, N>& m, const sm_vector<T, N>& w)
{
  sm_vector<T, N> erg;
  erg.assign_prod(m, w);
  return erg;
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_vector<T, N>& w, const sm_matrix_sym<T, N>& m)
{
  sm_vector<T, N> erg;
  erg.assign_prod(w, m);
  return erg;
}

//---------------------------------------------------------------------------------------
// Bilinear form
//---------------------------------------------------------------------------------------

template<class T, int N>
T sm_matrix_sym<T, N>::bilinear(const sm_vector<T, N>& v, const sm_vector<T, N>& w) const
{
  sm_vector<T, N> t;   //TODO: exploit symmetry
  t.assign_prod(*this, w);
  return v.dot(t);
}

template<class T, int N>
T sm_matrix_sym<T, N>::bilinear(const sm_vector<T, N>& v) const
{
  const me_t& m=*this;
  T erg = v(0) * v(0) * m(0,0);
  for (int i=1; i<N; i++) {
    T t = m(i,0)*v(0);
    for (int j=1; j<i; j++) {
      t += m(i,j)*v(j);
    }
    erg += v(i) * (v(i)*m(i,i) + 2*t);
  }
  return erg;
}
//---------------------------------------------------------------------------------------
// Determinant of symmetric 3-matrix
//---------------------------------------------------------------------------------------
template<class T>
T determinant(const sm_matrix_sym<T, 3> &m)
{
  T d = m(0,0) * m(1,1) * m(2,2)
        + 2 * m(0,1) * m(0,2) * m(1,2)
        - m(0,0) * m(1,2) * m(1,2)
        - m(1,1) * m(0,2) * m(0,2)
        - m(2,2) * m(0,1) * m(0,1);
  return d;
}

template<class T>
void invert_matrix(const sm_matrix_sym<T, 3>&m, sm_matrix_sym<T, 3>& erg, T& det)
{
  det = determinant(m);

  erg(0,0) = (-m(1,2)*m(1,2) + m(1,1)*m(2,2) );
  erg(0,1) = ( m(0,2)*m(1,2) - m(0,1)*m(2,2) );
  erg(1,1) = (-m(0,2)*m(0,2) + m(0,0)*m(2,2) );
  erg(0,2) = (-m(0,2)*m(1,1) + m(0,1)*m(1,2) );
  erg(1,2) = ( m(0,1)*m(0,2) - m(0,0)*m(1,2) );
  erg(2,2) = (-m(0,1)*m(0,1) + m(0,0)*m(1,1) );

  erg /= det;
}


//---------------------------------------------------------------------------------------
//  co/contra-variant vector
//---------------------------------------------------------------------------------------

template<class T, int N, bool UP> class sm_tensor1 {
  typedef sm_tensor1<T, N, UP> me_t;
  typedef sm_vector<T, N> vec_t;

  vec_t c;
  public:
  typedef T elem_t;
  enum {SIZE=N};

  sm_tensor1(){}
  sm_tensor1(zero_literal z) {c.zero();}

  T& operator()(int j) {return c(j);}
  const T& operator()(int j) const {return c(j);}
  vec_t& as_vector() {return c;}
  const vec_t& as_vector() const {return c;}

  void assign_sum(const me_t& a, const me_t& b) {c.assign_sum(a.c, b.c);}
  void assign_diff(const me_t& a, const me_t& b) {c.assign_diff(a.c, b.c);}
  void assign_minus(const me_t& a) {c.assign_minus(a.c);}
  void assign_prod(const me_t& a, const T& z) {c.assign_prod(a.c, z);}
  void assign_div(const me_t& a, const T& z) {c.assign_div(a.c, z);}
  template<bool UP1>
  void assign_prod(const sm_tensor1<T, N, !UP1>& v, const sm_tensor2_sym<T, N, UP1, UP>& m) {
    c.assign_prod(v.as_vector(), m.as_matrix());
  }
  template<bool UP1>
  void assign_prod(const sm_tensor2_sym<T, N, UP, UP1>& m, const sm_tensor1<T, N, !UP1>& v) {
    c.assign_prod(m.as_matrix(), v.as_vector());
  }

  void operator+=(const me_t &a) {assign_sum(*this, a);}
  void operator-=(const me_t &a) {assign_diff(*this, a);}
  void operator*=(const T& z) {assign_prod(*this, z);}
  void operator/=(const T& z) {assign_div(*this, z);}

  me_t operator+(const me_t &a) const {me_t e; e.assign_sum(*this, a); return e;}
  me_t operator-(const me_t &a) const {me_t e; e.assign_diff(*this, a); return e;}
  me_t operator*(const T& z) const {me_t e; e.assign_prod(*this, z); return e;}
  me_t operator/(const T& z) const {me_t e; e.assign_div(*this, z); return e;}
};

template<class T, int N, bool UP>
sm_tensor1<T, N, UP> operator-(sm_tensor1<T, N, UP> &v)
{
  sm_tensor1<T, N, UP> erg;
  erg.assign_minus(v);
  return erg;
}

//---------------------------------------------------------------------------------------
// Scalar * vector
//---------------------------------------------------------------------------------------

template<class T, int N, bool UP>
sm_tensor1<T, N, UP> inline operator*(const T& z, const sm_tensor1<T, N, UP>& a) {
  sm_tensor1<T, N, UP> erg;
  erg.assign_prod(a, z);
  return erg;
}

//---------------------------------------------------------------------------------------
// Vector-Vector contraction v^i w_i  resp. v_i w^i
//---------------------------------------------------------------------------------------

template<class T, int N, bool UP>
T contract(const sm_tensor1<T, N, UP>& v, const sm_tensor1<T, N, !UP>& w)
{
  return v.as_vector().dot(w.as_vector());
}

template<class T, int N, bool UP>
T operator*(const sm_tensor1<T, N, UP> &v, const sm_tensor1<T, N, !UP> &w)
{
  return contract(v, w);
}


//---------------------------------------------------------------------------------------
// Symmetric tensors m^ij , m_ij, m^i_j, m_i^j
//---------------------------------------------------------------------------------------


template<class T, int N, bool UP1, bool UP2> class sm_tensor2_sym {
  typedef sm_tensor2_sym<T, N, UP1, UP2> me_t;
  typedef sm_matrix_sym<T,N> mat_t;
  mat_t m;
  public:
  sm_tensor2_sym(){}
  sm_tensor2_sym(zero_literal z) {m.zero();}
  sm_tensor2_sym(one_literal z) {m.diag(1.0);}

  T& operator()(int i, int j) {return m(i,j);}
  const T& operator()(int i, int j) const {return m(i,j);}
  mat_t& as_matrix() {return m;}
  const mat_t& as_matrix() const {return m;}

  void assign_sum(const me_t& a, const me_t& b) {m.assign_sum(a.m, b.m);}
  void assign_diff(const me_t& a, const me_t& b) {m.assign_diff(a.m, b.m);}
  void assign_minus(const me_t& a) {m.assign_minus(a.m);}
  void assign_prod(const me_t& a, const T& z) {m.assign_prod(a.m, z);}
  void assign_div(const me_t& a, const T& z) {m.assign_div(a.m, z);}

  void operator+=(const me_t& a) {assign_sum(*this, a);}
  void operator-=(const me_t& a) {assign_diff(*this, a);}
  void operator*=(const T& z) {assign_prod(*this, z);}
  void operator/=(const T& z) {assign_div(*this, z);}

  me_t operator+(const me_t& a) const {me_t e; e.assign_sum(*this, a); return e;}
  me_t operator-(const me_t& a) const {me_t e; e.assign_diff(*this, a); return e;}
  me_t operator*(const T& z) const {me_t e; e.assign_prod(*this, z); return e;}
  me_t operator/(const T& z) const {me_t e; e.assign_div(*this, z); return e;}
  void diag(const T& d) {m.diag(d);}
  T contract(const sm_tensor1<T, N, !UP1>& v, const sm_tensor1<T, N, !UP2>& w) const {
    return m.bilinear(v.as_vector(), w.as_vector());
  }
  T quadratic(const sm_tensor1<T, N, !UP1>& v) const {
    return contract_quadratic(*this, v);
  }
};


template<class T, int N, bool UP1, bool UP2>
sm_tensor2_sym<T, N, UP1, UP2> operator-(sm_tensor2_sym<T, N, UP1, UP2> &v)
{
  sm_tensor2_sym<T, N, UP1, UP2> erg;
  erg.assign_minus(v);
  return erg;
}

//---------------------------------------------------------------------------------------
// Scalar * tensor
//---------------------------------------------------------------------------------------

template<class T, int N, bool UP1, bool UP2>
sm_tensor2_sym<T, N, UP1, UP2> operator*(const T& a, const sm_tensor2_sym<T, N, UP1, UP2> &m){
  sm_tensor2_sym<T, N, UP1, UP2> erg;
  erg.assign_prod(m, a);
  return erg;
}


//---------------------------------------------------------------------------------------
// Contraction m^ij w_j or m_ij w^j or m^i_j w^j
//---------------------------------------------------------------------------------------


template<class T, int N, bool UP1, bool UP2>
sm_tensor1<T, N, UP1> operator*(const sm_tensor2_sym<T, N, UP1, UP2>& m, const sm_tensor1<T, N, !UP2>& w)
{
  sm_tensor1<T, N, UP1> erg;
  erg.assign_prod(m, w);
  return erg;
}

template<class T, int N, bool UP1, bool UP2>
sm_tensor1<T, N, UP2> operator*(const sm_tensor1<T, N, !UP1>& w, const sm_tensor2_sym<T, N, UP1, UP2>& m)
{
  sm_tensor1<T, N, UP2> erg;
  erg.assign_prod(w, m);
  return erg;
}


//---------------------------------------------------------------------------------------
// Contraction v^i m_ij v^j
//---------------------------------------------------------------------------------------


template<class T, int N, bool UP>
T contract_quadratic(const sm_tensor2_sym<T, N, UP, UP>&m, const sm_tensor1<T, N, !UP> &v)
{
  return m.as_matrix().bilinear(v.as_vector());
}

//---------------------------------------------------------------------------------------
// Determinant
//---------------------------------------------------------------------------------------

template<class T, int N, bool UP1, bool UP2>
T determinant(const sm_tensor2_sym<T, N, UP1, UP2> &m)
{
  return determinant(m.as_matrix());
}

//---------------------------------------------------------------------------------------
// Metric
//---------------------------------------------------------------------------------------

template<class T, int N> struct sm_metric  {
  sm_tensor2_sym<T, N, false, false> lo;
  sm_tensor2_sym<T, N, true, true> up;
  T vol_elem, det;

  sm_metric(){}

  void raise(sm_tensor1<T, N, true>& erg, const sm_tensor1<T, N, false>& v) const {
    erg.assign_prod(up, v);
  }
  sm_tensor1<T, N, true> raise(const sm_tensor1<T, N, false>& v) const {
    sm_tensor1<T, N, true> erg;
    erg.assign_prod(up, v);
    return erg;
  }
  void lower(sm_tensor1<T, N, false>& erg, const sm_tensor1<T, N, true>& v) const {
    erg.assign_prod(lo, v);
  }
  sm_tensor1<T, N, false> lower(const sm_tensor1<T, N, true>& v) const {
    sm_tensor1<T, N, false> erg;
    erg.assign_prod(lo, v);
    return erg;
  }

  T contract(const sm_tensor1<T, N, true>& v, const sm_tensor1<T, N, true>& w) const {
    return lo.contract(v,w);
  }
  T contract(const sm_tensor1<T, N, false>& v, const sm_tensor1<T, N, false>& w) const {
    return up.contract(v,w);
  }
  T norm2(const sm_tensor1<T, N, true>& v) const {
    return lo.quadratic(v);
  }
  T norm2(const sm_tensor1<T, N, false>& v) const {
    return up.quadratic(v);
  }
  template<bool UP> T norm(const sm_tensor1<T, N, UP>& v) const {
    return std::sqrt(norm2(v));
  }

  void complete_from_lower() {
    invert_matrix(lo.as_matrix(), up.as_matrix(), det);
    vol_elem = std::sqrt(det);
  }
  void minkowski() {
    lo = ONE;
    up = ONE;
    vol_elem = 1.0;
    det      = 1.0;
  }
};

template<class T, bool UP1, bool UP2>
void collect(sm_tensor2_sym<T, 3, UP1, UP2>& m, const T& xx, const T& xy, const T& xz,
              const T& yy, const T& yz, const T& zz)
{
  m(0,0) = xx;
  m(0,1) = xy;
  m(0,2) = xz;
  m(1,1) = yy;
  m(1,2) = yz;
  m(2,2) = zz;
}

template<class T, bool UP1>
void collect(sm_tensor1<T, 3, UP1>& v, const T& x, const T& y, const T& z)
{
  v(0) = x;
  v(1) = y;
  v(2) = z;
}


}

#endif

