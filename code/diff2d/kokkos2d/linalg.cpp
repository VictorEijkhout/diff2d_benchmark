/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** linalg.cpp : bordered vector routines
 ****
 ****************************************************************/

#include <cmath>
#include <string>
using std::string;

#include <iostream>
using std::cout;
#include <format>
using std::format;

#include "linalg.hpp"

namespace linalg {
#include <algorithm>
  using std::for_each,std::transform;

  /*! Constructor with local allocation: 
    This creates a bordered array, storing the inner sizes
  */
  template< typename real >
  bordered_array<real>::bordered_array( size_t m,size_t n,int border )
    : _m(m),_n(n),border(border)
    , std::vector<real>( (m+2*border)*(n+2*border) ) {};

  //! Compute the 5-point Laplace stencil from an input array
  template< typename real >
  void bordered_array<real>::central_difference_from
      ( const linalg::bordered_array<real>& other,bool trace ) {
    for_each
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
        this->inner_range().begin(),this->inner_range().end(),
        [out = data2d(), in = other.data2d()] ( auto idx ) {
        auto [i,j] = idx;
        out[i,j] = 4*in[i,j] - in[i-1,j] - in[i+1,j] - in[i,j-1] - in[i,j+1]; }
        );
    log_flops(_m*_n*5); log_bytes( sizeof(real)*_m*_n*7 );
  };

  //! Copy the interior of another bordered array, but leave border alone
  template< typename real >
  void bordered_array<real>::copy_interior_from
      ( const linalg::bordered_array<real>& other ) {
    for_each
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
        this->inner_range().begin(),this->inner_range().end(),
        [out = data2d(),in = other.data2d()] ( auto idx ) {
        auto [i,j] = idx;
        out[i,j] = in[i,j]; }
        );
    log_flops(_m*_n*0); log_bytes( sizeof(real)*_m*_n*3 );
  };

  //! Scale the interior, leaving the border alone
  template< typename real >
  void bordered_array<real>::scale_interior
      ( const linalg::bordered_array<real>& other,real factor ) {
    const auto otherinner = other.inner_range();
    transform
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
	otherinner().begin(),otherinner().end(),
        this->inner_range().begin(),
        [out = data2d(),factor] ( auto idx ) {
        auto [i,j] = idx;
        out[i,j] *= factor; }
        );
    log_flops(_m*_n*1); log_bytes( sizeof(real)*_m*_n*2 );
  };

  //! Compute the L2 norm of the interior
  template< typename real >
  real bordered_array<real>::l2norm() {
    real sum_of_squares{0};
    for_each
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
        this->inner_range().begin(),this->inner_range().end(),
        [out = data2d(),&sum_of_squares] ( auto idx ) {
        auto [i,j] = idx;
        sum_of_squares += out[i,j] * out[i,j]; }
        );
    log_flops(_m*_n*3); log_bytes( sizeof(real)*_m*_n*1 );
    return std::sqrt(sum_of_squares);
  };

  //! Set the interior to a value
  template< typename real >
  void bordered_array<real>::set_value real value, bool trace ) {
    for_each
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
        this->inner_range().begin(),this->inner_range().end(),
        [this,out = data2d(),value,trace] ( auto idx ) {
        auto [i,j] = idx;
        out[i,j] = value;
      } );
    log_flops(_m*_n*0); log_bytes( sizeof(real)*_m*_n*2 );
  };

  template< typename real >
  void bordered_array<real>::set_bc(bool down, bool right, bool trace) {
    for_each
      ( 
#ifdef USE_TBB
        exec::par_unseq,
#endif
        domain().begin(),domain().end(),
        [this,out = data2d(),down,right] ( auto idx ) {
        auto [i,j] = idx;
        if ( down and i==_m-1 )
          out[i,j] += 1.;
        if ( right and j==_n-1 )
          out[i,j] += 1.;
      } );
  };

  //! Output the whole domain as cartesian array
  template< typename real >
  void bordered_array<real>::view( string caption ) {
    if (caption!="")
      cout << format("{}:\n",caption);
    for_each
      ( domain().begin(),domain().end(),
        [n=this->_n,array = data2d(),border=this->border] ( auto idx ) {
          auto [i,j] = idx;
          char c = ( j<n+2*border-1 ? ' ' : '\n' );
          cout << format("{:5.2}{}",array[i,j],c);
        }
        );
    cout << format("\n");
  };

  //! Copy internal data into a new vector
  template< typename real >
  std::vector<real> bordered_array<real>::internal_data() {
    std::vector<real> internal( _m * _n );
    size_t loc{0};
    std::for_each
      ( this->inner_range().begin(),this->inner_range().end(),
	[ out=internal.data(),in=data2d(),&loc ] ( auto idx ) {
	auto [i,j] = idx;
	out[loc++] = in[i,j];
      } );
    return internal;
  };
};

namespace linalg {
  template class bordered_array<float>;
  template class bordered_array<double>;
};

#if 0

double sum=0.;
Kokkos::parallel_reduce
( "ytAx product",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>( {0,0}, {m,n} ),
  KOKKOS_LAMBDA (int i,int j,double &partial ) {
  partial += yvec(i) * matrix(i,j) * xvec(j); },
        sum
  );

#endif
