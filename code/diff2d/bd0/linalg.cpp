/****************************************************************
 ****
 **** This file belongs with the course
 **** Introduction to Scientific Programming in C++/Fortran2003
 **** copyright 2019-2023 Victor Eijkhout eijkhout@tacc.utexas.edu
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

#include <fmt/format.h>
using fmt::print;

#include "../linalg.hpp"

#define IINDEX( i,j ) ((i)+border)*(_n+2*border) + (j)+border

namespace linalg {

  //! The constructor copies arguments and allocates the data
  template< typename real >
  bordered_array<real>::bordered_array( size_t m,size_t n,int border )
    : _m(m),_n(n),border(static_cast<size_t>(border))
    , std::vector<real>((m+2*border)*(n+2*border)) {}

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptomp
  template< typename real >
  void bordered_array<real>::central_difference_from
      ( linalg::bordered_array<real> other,bool trace ) {
    auto out = this->data();
    auto in = other.data();
    if (border==0) {
#   pragma omp parallel for collapse(clevel)
    for ( size_t i=1; i<_m-1; i++ ) {
      for ( size_t j=1; j<_n-1; j++ ) {
        out[ IINDEX(i,j) ] = 4*in[ IINDEX(i,j) ]
          - in[ IINDEX(i-1,j) ] - in[ IINDEX(i+1,j) ] - in[ IINDEX(i,j-1) ] - in[ IINDEX(i,j+1) ];
      }
    }
    } else {
#   pragma omp parallel for collapse(clevel)
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
  void bordered_array<real>::copy_interior_from( linalg::bordered_array<real> other ) {
    auto out = this->data();
    auto in = other.data();
#   pragma omp parallel for collapse(clevel)
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
  void bordered_array<real>::scale_interior( real factor ) {
    auto out = this->data();
#   pragma omp parallel for collapse(clevel)
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ )
        out[ IINDEX(i,j) ] *= factor;
    log_flops(_m*_n*1); log_bytes( sizeof(real)*_m*_n*2 );
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormomp
  template< typename real >
  real bordered_array<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
#   pragma omp parallel for collapse(clevel) reduction(+:sum_of_squares)
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
  void bordered_array<real>::set( real value, bool trace ) {
    auto out = this->data();
#   pragma omp parallel for collapse(clevel)
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ ) {
        auto ij = IINDEX(i,j);
        out[ ij ] = value;
      }
    log_flops(_m*_n*0); log_bytes( sizeof(real)*_m*_n*2 );
  };

  template< typename real >
  void bordered_array<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
#   pragma omp parallel for collapse(clevel)
    for ( size_t i=0; i<_m; i++ )
      for ( size_t j=0; j<_n; j++ )
        if ( i==_m-1 or j==_n-1 )
          out[ IINDEX(i,j) ] = 1.;
  };

  template< typename real >
  void bordered_array<real>::view( string caption ) {
    if (caption!="")
      print("{}:\n",caption);
    auto out = this->data();
    for ( size_t i=0; i<_m+2*border; i++ ) {
      for ( size_t j=0; j<_n+2*border; j++ ) {
        char c = ( j<_n+2*border-1 ? ' ' : '\n' );
        print("{:5.2}{}",out[ oindex(i,j) ],c);
      }
    }
  };
};

namespace linalg {
  template class bordered_array<float>;
  template class bordered_array<double>;
};
