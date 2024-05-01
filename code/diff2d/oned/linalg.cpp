/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** linalg.cpp : bordered vector routines for old-style omp
 ****
 ****************************************************************/

#define USE_OMP

#include <cmath>
#include <algorithm>
using std::for_each;
#include <string>
using std::string;

// #include <fmt/format.h>
// using fmt::print;
#include <iostream>
using std::cout;
#include <format>
using std::format;

#include "../linalg.hpp"

#define IINDEX( i,j ) ((i)+border)*_n2b + (j)+border

namespace linalg {

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptomp
  template< typename real >
  void bordered_array_1d<real>::central_difference_from
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    if (border==0) {
#   pragma omp parallel for 
    for ( size_t i=1; i<_m-1; i++ ) {
      for ( size_t j=1; j<_n-1; j++ ) {
        out[ IINDEX(i,j) ] = 4*in[ IINDEX(i,j) ]
          - in[ IINDEX(i-1,j) ] - in[ IINDEX(i+1,j) ] - in[ IINDEX(i,j-1) ] - in[ IINDEX(i,j+1) ];
      }
    }
    } else {
#   pragma omp parallel for 
    for ( size_t i=0; i<_m; i++ ) {
      for ( size_t j=0; j<_n; j++ ) {
        out[ IINDEX(i,j) ] = 4*in[ IINDEX(i,j) ]
          - in[ IINDEX(i-1,j) ] - in[ IINDEX(i+1,j) ] - in[ IINDEX(i,j-1) ] - in[ IINDEX(i,j+1) ];
      }
    }
    }
    log_flops(_m*_n*5); log_bytes( sizeof(real)*_m*_n*7 );
  };
  //codesnippet end

  //! Copy the interior of another bordered array, but leave border alone
  //codesnippet d2dcopyomp
  template< typename real >
  void bordered_array_1d<real>::copy_interior_from
      ( const linalg::bordered_array_base<real>& _other ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
#   pragma omp parallel for 
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ ) {
        out[ IINDEX(i,j) ] = in[ IINDEX(i,j) ];
      }
    log_flops(_m*_n*0); log_bytes( sizeof(real)*_m*_n*3 );
  };
  //codesnippet end

  //! Scale the interior, leaving the border alone
  //codesnippet d2dscaleomp
  template< typename real >
  void bordered_array_1d<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
#   pragma omp parallel for 
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ )
        out[ IINDEX(i,j) ] = in[ IINDEX(i,j) ] * factor;
    log_flops(_m*_n*1); log_bytes( sizeof(real)*_m*_n*2 );
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormomp
  template< typename real >
  real bordered_array_1d<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
#   pragma omp parallel for reduction(+:sum_of_squares)
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ ) {
        auto v = out[ IINDEX(i,j) ];
        sum_of_squares += v*v;
      }
    log_flops(_m*_n*3); log_bytes( sizeof(real)*_m*_n*1 );
    return std::sqrt(sum_of_squares);
  };
  //codesnippet end

  //! Set the interior to a value
  template< typename real >
  void bordered_array_1d<real>::set( real value, bool trace ) {
    auto out = this->data();
#   pragma omp parallel for 
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ ) {
        auto ij = IINDEX(i,j);
        out[ ij ] = value;
      }
    log_flops(_m*_n*0); log_bytes( sizeof(real)*_m*_n*2 );
  };

  template< typename real >
  void bordered_array_1d<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
#   pragma omp parallel for 
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ )
        if ( i==_m-1 or j==_n-1 )
          out[ IINDEX(i,j) ] = 1.;
  };

  template< typename real >
  void bordered_array_1d<real>::view( string caption ) {
    if (caption!="")
      cout << format("{}:\n",caption);
    auto out = this->data();
    for ( size_t i=0; i<_m+2*border; i++ ) {
      for ( size_t j=0; j<_n+2*border; j++ ) {
        char c = ( j<_n+2*border-1 ? ' ' : '\n' );
        cout << format("{:5.2}{}",out[ oindex(i,j) ],c);
      }
    }
  };

  //! Copy internal data into a new vector
  template< typename real >
  std::vector<real> bordered_array_1d<real>::internal_data() {
    std::vector<real> internal( _m * _n );
    size_t loc{0};
    std::for_each
      ( this->inner().begin(),this->inner().end(),
	[ out=internal.data(),in=data2d(),&loc ] ( auto idx ) {
	auto [i,j] = idx;
	out[loc++] = in[i,j];
      } );
    return internal;
  };
};

namespace linalg {
  template class bordered_array_1d<float>;
  template class bordered_array_1d<double>;
};
