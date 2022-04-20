#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

#include <string>
#include <cmath>
#include <complex>
#include <iostream>
using namespace std;

typedef double DP;

template <class T>
class NRVec {
 private:
  int nn;// size of array. upper index is nn-1
  T *v;
 public:
  NRVec();
  explicit NRVec(int n);// Zero-based array
  NRVec(const T &a, int n);//initialize to constant value
  NRVec(const T *a, int n);// Initialize to array
  NRVec(const NRVec &rhs);// Copy constructor
  NRVec & operator=(const NRVec &rhs);//assignment
  NRVec & operator=(const T &a);//assign a to every element
  inline T & operator[](const int i);//i'th element
  inline const T & operator[](const int i) const;
  inline int size() const;
  ~NRVec();
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) {}

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(new T[n])
{
  for(int i=0; i<n; i++)
    v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
  for(int i=0; i<n; i++)
    v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
  for(int i=0; i<nn; i++)
    v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
     // postcondition: normal assignment via copying has been performed;
     //if vector and rhs were different sizes, vector
     //has been resized to match the size of rhs
{
  if (this != &rhs)
    {
      if (nn != rhs.nn) {
	if (v != 0) delete [] (v);
	nn=rhs.nn;
	v= new T[nn];
      }
      for (int i=0; i<nn; i++)
	v[i]=rhs[i];
    }
  return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)//assign a to every element
{
  for (int i=0; i<nn; i++)
    v[i]=a;
  return *this;
}

template <class T>
inline T & NRVec<T>::operator[](const int i)//subscripting
{
  return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const//subscripting
{
  return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
  return nn;
}

template <class T>
NRVec<T>::~NRVec()
{
  if (v != 0)
    delete[] (v);
}



#endif _NR_UTIL_H_
