/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** seq.cpp : one-dimensional indexing for diff2d
 ****
 ****************************************************************/

#include <cmath>
#include <format>
#include <iostream>
#include <string>

#include "seq.hpp"

//codesnippet d2dseqindex
#define IINDEX( i,j,m,n,b ) ((i)+b)*(n+2*b) + (j)+b
//codesnippet end

namespace sparsealg {

  template< typename real >
  bordered_array_seq<real>::bordered_array_seq( idxint m_inner,idxint n_inner,int border )
    : bordered_array_base<real>(m_inner,n_inner,border) {
    auto out = this->data();
    const auto [m,n,b] = this->outer_sizes();
    //std::cout << std::format("zero with border: {} x {}\n",m,n);
    for ( idxint i=0; i<m; ++i )
      for ( idxint j=0; j<n; ++j )
        out[ index2d(i,j,m,n,0) ] = static_cast<real>(0);
  };

  //! Compute the 5-point Laplace stencil from an input array
  template< typename real >
  void bordered_array_seq<real>::central_difference_from
      ( const sparsealg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const sparsealg::bordered_array_seq<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
  //codesnippet d2d5ptseq
    for ( idxint i=0; i<m; i++ ) {
      for ( idxint j=0; j<n; j++ ) {
        out[ index2d(i,j,m,n,b) ] = 4*in[ index2d(i,j,m,n,b) ]
          - in[ index2d(i-1,j,m,n,b) ] - in[ index2d(i+1,j,m,n,b) ]
          - in[ index2d(i,j-1,m,n,b) ] - in[ index2d(i,j+1,m,n,b) ];
      }
    }
    //codesnippet end
    log_flops(m*n*5); log_bytes( sizeof(real)*m*n*7 );
  };

  //! Scale the interior, leaving the border alone
  //codesnippet d2dupcast
  template< typename real >
  void bordered_array_seq<real>::scale_interior
      ( const sparsealg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other =
      dynamic_cast<const sparsealg::bordered_array_seq<real>&>(_other);
  //codesnippet end
    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    //codesnippet d2dscaleseq
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ )
        out[ index2d(i,j,m,n,b) ] = in[ index2d(i,j,m,n,b) ] * factor;
    //codesnippet end
    log_flops(m*n*1); log_bytes( sizeof(real)*m*n*2 );
  };

  //! Compute the L2 norm of the interior
  template< typename real >
  real bordered_array_seq<real>::l2norm() {
    real sum_of_squares{0};
    auto in = this->data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
  //codesnippet d2dnormseq
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ ) {
        auto v = in[ index2d(i,j,m,n,b) ];
        sum_of_squares += v*v;
      }
    log_flops(m*n*3); log_bytes( sizeof(real)*m*n*1 ); //snippetskip
    return std::sqrt(sum_of_squares);
  //codesnippet end
  };

  //! Set the interior to a value
  template< typename real >
  void bordered_array_seq<real>::set_value( real value, bool trace ) {
    auto out = this->data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ ) {
        auto ij = index2d(i,j,m,n,b);
        out[ ij ] = value;
      }
    log_flops(m*n*0); log_bytes( sizeof(real)*m*n*2 );
  };

  template< typename real >
  void bordered_array_seq<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ )
        if ( i==m-1 or j==n-1 )
          out[ index2d(i,j,m,n,b) ] = 1.;
  };

};

namespace sparsealg {
  template class bordered_array_seq<float>;
  template class bordered_array_seq<double>;
};
