/**********************************************************************
 * flux: framework for learning unstructured meshing
 * Copyright (c) 2022 Philip Caplan. All rights reserved.
 * Licensed under the MIT License (https://mit-license.org/)
 **********************************************************************/
#ifndef FLUX_VEC_H_
#define FLUX_VEC_H_

#include "error.h"
#include "mat.h"
#include "result_of.h"

#include <iostream>
#include <vector>

namespace flux {

/**
 * \brief Represents a m-dimensional vector whose entries are allocated dynamically.
 *        Inherits from matd, thus it is a m x 1 (m rows, 1 column) matrix.
 *        Templated by type T so the vector could store any type of entry,
 *        provided arithmetic operators are defined for the entry type.
 */
template<typename T>
class vecd : public matd<T> {
public:

  /**
   * \brief Constructs an m-dimensional vector.
   *
   * \param[in] m - number of components in the vector
   */
  vecd(int m) :
    matd<T>(m,1)
  {}

  /**
   * \brief Constructs an m-dimensional vector from some array of values.
   *
   * \param[in] m - number of components in the vector
   * \param[in] x - pointer to array of values to set into this vector
   */
  vecd( int m , const T* x ) :
    matd<T>(m,1) {
    for (int i = 0; i < m; i++)
      data_[i] = x[i];
  }

  /**
   * \brief Constructs an m-dimensional vector from some array of values.
   *
   * \param[in] x - std::vector of values to set into this vector (m = x.size())
   */
  vecd( const std::vector<T>& x ) :
    matd<T>(x.size(),1) {
    for (int i = 0; i < m_; i++)
      data_[i] = x[i];
  }

  /**
   * \brief Sets the values stored in some vector into this vector.
   *
   * \param[in] x - vector to copy
   */
  void set( const vecd<T>& x ) {
    flux_assert( x.m() == m_ );
    for (int i = 0; i < m_; i++)
      data_[i] = x(i);
  }

  /**
   * \brief Read/write acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return reference to entry at i'th index
   */
  T& operator() (int i) {
    flux_assert_msg( i < m_ , "attempt to access i = %d but m = %d" , i , m_ );
    return data_[i];
  }

  /**
   * \brief Read-only acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return const reference to entry at i'th index
   */
  const T& operator() (int i) const {
    flux_assert_msg( i < m_ , "attempt to access i = %d but m = %d" , i , m_ );
    return data_[i];
  }

  /**
   * \brief Read/write acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return reference to entry at i'th index
   */
  T& operator[] (int i) {
    flux_assert( i < m_ );
    return data_[i];
  }

  /**
   * \brief Read-only acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return const reference to entry at i'th index
   */
  const T& operator[] (int i) const {
    flux_assert( i < m_ );
    return data_[i];
  }

  /**
   * \brief Set all entries in this vector to zero.
   */
  void zero() {
    for (int i = 0; i < m_; i++)
      data_[i] = 0;
  }

  /**
   * \brief Return the number of components in this vector.
   */
  int m() const { return m_; }

  /**
   * \brief Print out the components of the vector
   *
   * \param[in] v0 (optional) - prefix for printing vector entries
   */
  void print( const std::string& v0 = std::string() ) const {
    std::string v = v0.empty() ? "" : v0;
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (int i = 0; i < m_; i++)
      std::cout << v << "(" + std::to_string(i) + "): " << (*this)(i) << std::endl;
  }

  /**
   * \brief Read/write access to the data vector storing the components of this vector.
   *
   * \return pointer to first element in the data array
   */
  T* data() { return data_.data(); }

  /**
   * \brief Read-only access to the data vector storing the components of this vector.
   *
   * \return pointer to first element in the data array
   */
  const T* data() const { return data_.data(); }

private:
  using matd<T>::m_;    // number of components in the vector
  using matd<T>::data_; // array of vector components
};

/**
 * \brief Represents a m-dimensional vector whose entries are allocated statically.
 *        Inherits from mats, thus it is a m x 1 (m rows, 1 column) matrix.
 *        Templated by type T so the vector could store any type of entry,
 *        provided arithmetic operators are defined for the entry type.
 */
template<int _M,typename T>
class vecs : public mats<_M,1,T> {

public:
  static const int M = _M;

private:
  using mats<M,1,T>::data_;

public:
  /**
   * \brief Constructs an M-dimensional vector
   */
  vecs() {}

  /**
   * \brief Constructs an m-dimensional vector from some array of values,
   *        checking that _m is <= M
   *
   * \param[in] data - pointer to array of values to set into this vector
   * \param[in] _m - number of components to save into the vector (should be <= M)
   */
  vecs( const T* data , int _m = -1 ) {
    flux_assert( _m <= M );
    if (_m < 0) _m = M;
    for (int i = 0; i < _m; i++)
      data_[i] = data[i];
  }

  /**
   * \brief Constructs an M-dimensional vector and sets all entries to the integer a.
   *
   * \param[in] a - integer to set all components to
   */
  vecs( int a ) {
    for (int i = 0; i < M; i++)
      data_[i] = a;
  }

  /**
   * \brief Copies the components from some other vector b into this vector.
   *
   * \param[in] b - vector to copy of a potentially different type
   */
  template<typename S>
  vecs(const vecs<M,S>& b) {
    for (int i = 0; i < M; i++)
      (*this)(i) = b(i);
  }

  /**
   * \brief Copies the components from some other M x 1 matrix into this vector.
   *
   * \param[in] A - vector (M x 1 matrix) to copy
   */
  vecs( const mats<M,1,T>& A ) {
    for (int i = 0; i < M; i++)
      data_[i] = A(i,0);
  }

  /**
   * \brief Assigns the components from some other vector b into this vector.
   *        Useful for writing expressions such as vecs<3,double> u = v.
   *
   * \param[in] b - vector to copy of a potentially different type
   */
  template<typename S>
  vecs<M,typename result_of<S,T>::type>& operator= (const vecs<M,S>& b) {
    for (int i = 0; i < M; i++)
      (*this)(i) = b(i);
    return *this;
  }

  /**
   * \brief Assigns the components from some initializer list into this vector.
   *        Useful for writing thing like vecs<3,double> x = {1,2,3}.
   *
   * \param[in] v - initializer list to copy of a potentially different type
   */
  template<typename S>
  vecs<M,typename result_of<S,T>::type>& operator= (const std::initializer_list<S>& v) {
    flux_assert( v.size() == M );
    int i = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      (*this)(i++) = *it;
    return *this;
  }

  /**
   * \brief Assigns the components from some initializer list into this vector.
   *        Useful for writing thing like vecs<3,double> x = {1,2,3}.
   *
   * \param[in] v - initializer list to copy of a potentially different type
   */
  template<typename R>
  vecs( const std::initializer_list<R>& v ) {
    operator=( v );
  }

  /**
   * \brief Assigns the components from some initializer list into this vector.
   *        Useful for writing thing like vecs<3,double> x = {1,2,3}.
   *
   * \param[in] v - initializer list to copy of a potentially different type
   */
  vecs<M,T>& operator= (const std::initializer_list<T>& v ) {
    flux_assert( v.size() == M );
    int i = 0;
    for (auto it = v.begin(); it != v.end(); ++it)
      (*this)(i++) = *it;
    return *this;
  }

  /**
   * \brief Sets all components of this vector to a scalar value.
   *
   * \param[in] a - scalar to set all entries to
   */
  vecs<M,T>& operator= (int a) {
    for (int i = 0; i < M; i++)
      data_[i] = a;
    return *this;
  }

  /**
   * \brief Read/write acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return reference to entry at i'th index
   */
  T& operator() (int i) {
    flux_assert( i < M );
    return data_[i];
  }

  /**
   * \brief Read-only acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return const reference to entry at i'th index
   */
  const T& operator() (int i) const {
    flux_assert( i < M );
    return data_[i];
  }

  /**
   * \brief Read/write acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return reference to entry at i'th index
   */
  T& operator[] (int i) {
    flux_assert( i < M );
    return data_[i];
  }

  /**
   * \brief Read-only acces to i'th component.
   *
   * \param[in] i - index of component to access
   *
   * \return const reference to entry at i'th index
   */
  const T& operator[] (int i) const {
    flux_assert( i < M );
    return data_[i];
  }

  /**
   * \brief Print out the components of the vector
   *
   * \param[in] v0 (optional) - prefix for printing vector entries
   */
  void print( const std::string& v0 = std::string() ) const {
    std::string v = v0.empty() ? "" : v0;
    printf("%s:\n",__PRETTY_FUNCTION__);
    for (int i = 0; i < M; i++)
      std::cout << v << "(" + std::to_string(i) + "): " << (*this)(i) << std::endl;
  }
};

/**
 * \brief Computes the dot product u.v between two static vectors.
 *
 * \param[in] u - vector (vecs)
 * \param[in] v - vector (vecs)
 *
 * \return dot product u.v
 */
template<typename R,typename S,int M> typename result_of<R,S>::type dot( const vecs<M,R>& u , const vecs<M,S>& v );

/**
 * \brief Computes the dot product u.v between two dynamic vectors.
 *
 * \param[in] u - vector (vecd)
 * \param[in] v - vector (vecd)
 *
 * \return dot product u.v
 */
template<typename R,typename S> typename result_of<R,S>::type dot( const vecd<R>& u , const vecd<S>& v );

/**
 * \brief Normalizes a vector so it has a unit length.
 *
 * \param[in,out] u - vector (static) to normalize
 */
template<typename T, int M> void normalize( vecs<M,T>& u );

/**
 * \brief Computes the magnitude (length) of a static vector
 *
 * \param[in] u - vector to compute the length of
 *
 * \return length of u: || u || = sqrt( u^T u )
 */
template<typename T, int M> T magnitude( const vecs<M,T>& u );

/**
 * \brief Computes the magnitude (length) of a dynamic vector
 *
 * \param[in] u - vector to compute the length of
 *
 * \return length of u: || u || = sqrt( u^T u )
 */
template<typename T> T magnitude( const vecd<T>& u );

/**
 * \brief Computes the cross-product of two 3d vectors: w = u x v
 *
 * \param[in] u - 3d vector (left operand)
 * \param[in] v - 3d vector (right operand)
 *
 * \return 3d vector u x v
 */
template<typename T> vecs<3,T> cross( const vecs<3,T>& u , const vecs<3,T>& v );

/**
 * \brief Computes the (dynamic) vector subtraction x - y
 */
template<typename R,typename S> vecd< typename result_of<R,S>::type > operator-( const vecd<R>& x , const vecd<S>& y );

/**
 * \brief Computes the (dynamic) vector addition x + y
 */
template<typename R,typename S> vecd< typename result_of<R,S>::type > operator+( const vecd<R>& x , const vecd<S>& y );

/**
 * \brief Computes the scalar-vector multiplication a * x
 */
template<typename R,typename S> vecd< typename result_of<R,S>::type > operator*( const R& a , const vecd<S>& y );

/**
 * \brief Computes the vector-scalar multiplication x * a
 */
template<typename R,typename S> vecd< typename result_of<R,S>::type > operator*( const vecd<R>& x , const S& a );

/**
 * \brief Computes the vector-scalar division x / a
 */
template<typename R,typename S> vecd< typename result_of<R,S>::type > operator/( const vecd<R>& x , const S& a );

/**
 * \brief Specialized definition for 2d static double vectors.
 */
typedef vecs<2,double> vec2d;

/**
 * \brief Specialized definition for 3d static double vectors.
 */
typedef vecs<3,double> vec3d;

/**
* \brief Specialized definition for 3d static float vectors.
 */
typedef vecs<3,float>  vec3f;

/**
* \brief Specialized definition for 4d static double vectors.
 */
typedef vecs<4,double> vec4d;

} // flux

#endif
