/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** seq.cpp : one-dimensional indexing for diff2d
 ****
 ****************************************************************/

#include <cmath>
#include <format>
#include <iostream>
#include <string>

#include "seq.hpp"

//codesnippet d2d5ptompindex
#define IINDEX( i,j,b,n2b ) ((i)+b)*n2b + (j)+b
//codesnippet end

namespace linalg {

  template< typename real >
  bordered_array_seq<real>::bordered_array_seq( size_t m,size_t n,int border )
    : bordered_array_base<real>(m,n,border) {
    auto out = this->data();
    auto b = this->border(), n2b = this->n2b();
    for ( auto ij=0; ij<(m+2*b)*(n+2*b); ++ij )
      out[ij] = static_cast<real>(0);
  };

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptomp
  template< typename real >
  void bordered_array_seq<real>::central_difference_from
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_seq<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    std::cout << std::format("seq m={} n={} b={} n2b={}\n",m,n,b,n2b);
    for ( size_t i=0; i<m; i++ ) {
      for ( size_t j=0; j<n; j++ ) {
        out[ IINDEX(i,j,b,n2b) ] = 4*in[ IINDEX(i,j,b,n2b) ]
          - in[ IINDEX(i-1,j,b,n2b) ] - in[ IINDEX(i+1,j,b,n2b) ] - in[ IINDEX(i,j-1,b,n2b) ] - in[ IINDEX(i,j+1,b,n2b) ];
      }
    }
  //codesnippet end
    log_flops(m*n*5); log_bytes( sizeof(real)*m*n*7 );
  };

  //! Scale the interior, leaving the border alone
  //codesnippet d2dscaleomp
  template< typename real >
  void bordered_array_seq<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_seq<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    for ( size_t i=0; i<m; i++ )
      for ( size_t j=0; j<n; j++ )
        out[ IINDEX(i,j,b,n2b) ] = in[ IINDEX(i,j,b,n2b) ] * factor;
    log_flops(m*n*1); log_bytes( sizeof(real)*m*n*2 );
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormomp
  template< typename real >
  real bordered_array_seq<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    for ( size_t i=0; i<m; i++ )
      for ( size_t j=0; j<n; j++ ) {
        auto v = out[ IINDEX(i,j,b,n2b) ];
        sum_of_squares += v*v;
      }
    log_flops(m*n*3); log_bytes( sizeof(real)*m*n*1 );
    return std::sqrt(sum_of_squares);
  };
  //codesnippet end

  //! Set the interior to a value
  template< typename real >
  void bordered_array_seq<real>::set( real value, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    for ( size_t i=0; i<m; i++ )
      for ( size_t j=0; j<n; j++ ) {
        auto ij = IINDEX(i,j,b,n2b);
        out[ ij ] = value;
      }
    log_flops(m*n*0); log_bytes( sizeof(real)*m*n*2 );
  };

  template< typename real >
  void bordered_array_seq<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    for ( size_t i=0; i<m; i++ )
      for ( size_t j=0; j<n; j++ )
        if ( i==m-1 or j==n-1 )
          out[ IINDEX(i,j,b,n2b) ] = 1.;
  };

};

namespace linalg {
  template class bordered_array_seq<float>;
  template class bordered_array_seq<double>;
};
