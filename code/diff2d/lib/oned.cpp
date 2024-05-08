/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** oned.cpp : one-dimensional indexing for diff2d
 ****
 ****************************************************************/

#include <cmath>
#include <format>
#include <iostream>
#include <string>

#include "omp.h"
#include "oned.hpp"

//codesnippet d2donedindex
#define IINDEX( i,j ) ((i)+border)*n2b + (j)+border
//codesnippet end

namespace linalg {

  template< typename real >
  bordered_array_1d<real>::bordered_array_1d( int64_t m,int64_t n,int border )
    : bordered_array_base<real>(m,n,border) {
    auto out = this->data();
    auto b = this->border();
    auto n2b = this->n2b();
    #pragma omp parallel for
    for ( auto ij=0; ij<(m+2*b)*(n+2*b); ++ij )
      out[ij] = static_cast<real>(0);
  };

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptoned
  template< typename real >
  void bordered_array_1d<real>::central_difference_from
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    #pragma omp parallel for
    for ( int64_t i=0; i<m; i++ ) {
      for ( int64_t j=0; j<n; j++ ) {
        out[ IINDEX(i,j) ] = 4*in[ IINDEX(i,j) ]
          - in[ IINDEX(i-1,j) ] - in[ IINDEX(i+1,j) ]
          - in[ IINDEX(i,j-1) ] - in[ IINDEX(i,j+1) ];
      }
    }
  //codesnippet end
    log_flops(m*n*5); log_bytes( sizeof(real)*m*n*3 );
  };

  //! Scale the interior, leaving the border alone
  //codesnippet d2dscaleoned
  template< typename real >
  void bordered_array_1d<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    #pragma omp parallel for 
    for ( int64_t i=0; i<m; i++ )
      for ( int64_t j=0; j<n; j++ )
        out[ IINDEX(i,j) ] = in[ IINDEX(i,j) ] * factor;
    log_flops(m*n*1); log_bytes( sizeof(real)*m*n*3 );
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormomp
  template< typename real >
  real bordered_array_1d<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    #pragma omp parallel for reduction(+:sum_of_squares)
    for ( int64_t i=0; i<m; i++ )
      for ( int64_t j=0; j<n; j++ ) {
        auto v = out[ IINDEX(i,j) ];
        sum_of_squares += v*v;
      }
    log_flops(m*n*3); log_bytes( sizeof(real)*m*n*1 );
    return std::sqrt(sum_of_squares);
  };
  //codesnippet end

  //! Set the interior to a value
  template< typename real >
  void bordered_array_1d<real>::set( real value, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    #pragma omp parallel for 
    for ( int64_t i=0; i<m; i++ )
      for ( int64_t j=0; j<n; j++ ) {
        auto ij = IINDEX(i,j);
        out[ ij ] = value;
      }
    log_flops(m*n*0); log_bytes( sizeof(real)*m*n*2 );
  };

  template< typename real >
  void bordered_array_1d<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    #pragma omp parallel for 
    for ( int64_t i=0; i<m; i++ )
      for ( int64_t j=0; j<n; j++ )
        if ( i==m-1 or j==n-1 )
          out[ IINDEX(i,j) ] = 1.;
  };

  template< typename real >
  void bordered_array_1d<real>::view( std::string caption ) {
    if (caption!="")
      std::cout << format("{}:\n",caption);
    auto out = this->data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    for ( int64_t i=0; i<m+2*border; i++ ) {
      for ( int64_t j=0; j<n+2*border; j++ ) {
        char c = ( j<n+2*border-1 ? ' ' : '\n' );
        std::cout << std::format("{:5.2}{}",out[ this->oindex(i,j) ],c);
      }
    }
  };

};

namespace linalg {
  template class bordered_array_1d<float>;
  template class bordered_array_1d<double>;
};
