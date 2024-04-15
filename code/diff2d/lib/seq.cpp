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

//codesnippet d2dindexscale
#define IINDEX( i,j,b,n2b ) ((i)+b)*n2b + (j)+b
//codesnippet end

namespace linalg {

  template< typename real >
  bordered_array_seq<real>::bordered_array_seq( int64_t m_inner,int64_t n_inner,int border )
    : bordered_array_base<real>(m_inner,n_inner,border) {
    auto out = this->data();
    const auto [m,n,b] = this->outer_sizes();
    //std::cout << std::format("zero with border: {} x {}\n",m,n);
    for ( int64_t i=0; i<m; ++i )
      for ( int64_t j=0; j<n; ++j )
	out[ IINDEX(i,j,0,n) ] = static_cast<real>(0);
  };

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptseq
  template< typename real >
  void bordered_array_seq<real>::central_difference_from
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_seq<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    for ( int64_t i=0; i<m; i++ ) {
      for ( int64_t j=0; j<n; j++ ) {
        out[ IINDEX(i,j,b,n2b) ] = 4*in[ IINDEX(i,j,b,n2b) ]
          - in[ IINDEX(i-1,j,b,n2b) ] - in[ IINDEX(i+1,j,b,n2b) ]
	  - in[ IINDEX(i,j-1,b,n2b) ] - in[ IINDEX(i,j+1,b,n2b) ];
      }
    }
  //codesnippet end
    log_flops(m*n*5); log_bytes( sizeof(real)*m*n*7 );
  };

  //! Scale the interior, leaving the border alone
  template< typename real >
  void bordered_array_seq<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const linalg::bordered_array_seq<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    //codesnippet d2dscaledseq
    for ( size_t i=0; i<m; i++ )
      for ( size_t j=0; j<n; j++ )
        out[ IINDEX(i,j,b,n2b) ] = in[ IINDEX(i,j,b,n2b) ] * factor;
    //codesnippet end
    log_flops(m*n*1); log_bytes( sizeof(real)*m*n*2 );
  };

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormseq
  template< typename real >
  real bordered_array_seq<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
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
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
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
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
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
