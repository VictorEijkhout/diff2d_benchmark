/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** clps.cpp : collapsed OpenMP
 ****
 ****************************************************************/

#include <cmath>
#include <format>
#include <iostream>
#include <string>

#include "omp.h"
#include "clps.hpp"

#define IINDEX( i,j,m,n,b ) ((i)+b)*(n+2*b) + (j)+b

namespace sparsealg {

  template< typename real >
  bordered_array_1d<real>::bordered_array_1d( idxint m,idxint n,int border )
    : bordered_array_base<real>(m,n,border) {
    auto out = this->data();
    auto b = this->border();
    auto n2b = this->n2b();
    #pragma omp parallel for
    for ( auto ij=0; ij<(m+2*b)*(n+2*b); ++ij )
      out[ij] = static_cast<real>(0);
  };

  //! Compute the 5-point Laplace stencil from an input array
  template< typename real >
  void bordered_array_1d<real>::central_difference_from
      ( const sparsealg::bordered_array_base<real>& _other,bool trace ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const sparsealg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    #pragma omp parallel for collapse(2)
    for ( idxint i=0; i<m; i++ ) {
      for ( idxint j=0; j<n; j++ ) {
        out[ IINDEX(i,j,m,n,b) ] = 4*in[ IINDEX(i,j,m,n,b) ]
          - in[ IINDEX(i-1,j,m,n,b) ] - in[ IINDEX(i+1,j,m,n,b) ]
          - in[ IINDEX(i,j-1,m,n,b) ] - in[ IINDEX(i,j+1,m,n,b) ];
      }
    }
    log_flops(m*n*5); log_bytes( sizeof(real)*m*n*3 );
  };

  //! Scale the interior, leaving the border alone
  template< typename real >
  void bordered_array_1d<real>::scale_interior
      ( const sparsealg::bordered_array_base<real>& _other, real factor ) {
    // upcast base to derived type
    const auto& other = dynamic_cast<const sparsealg::bordered_array_1d<real>&>(_other);
    auto out = this->data();
    auto in = other.data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    #pragma omp parallel for collapse(2)
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ )
        out[ IINDEX(i,j,m,n,b) ] = in[ IINDEX(i,j,m,n,b) ] * factor;
    log_flops(m*n*1); log_bytes( sizeof(real)*m*n*3 );
  };

  //! Compute the L2 norm of the interior
  template< typename real >
  real bordered_array_1d<real>::l2norm() {
    real sum_of_squares{0};
    auto out = this->data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    #pragma omp parallel for collapse(2) reduction(+:sum_of_squares)
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ ) {
        auto v = out[ IINDEX(i,j,m,n,b) ];
        sum_of_squares += v*v;
      }
    log_flops(m*n*3); log_bytes( sizeof(real)*m*n*1 );
    return std::sqrt(sum_of_squares);
  };

  //! Set the interior to a value
  template< typename real >
  void bordered_array_1d<real>::set_value( real value, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    #pragma omp parallel for collapse(2)
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ ) {
        auto ij = IINDEX(i,j,m,n,b);
        out[ ij ] = value;
      }
    log_flops(m*n*0); log_bytes( sizeof(real)*m*n*2 );
  };

  template< typename real >
  void bordered_array_1d<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    #pragma omp parallel for 
    for ( idxint i=0; i<m; i++ )
      for ( idxint j=0; j<n; j++ )
        if ( i==m-1 or j==n-1 )
          out[ IINDEX(i,j,m,n,b) ] = 1.;
  };

  template< typename real >
  void bordered_array_1d<real>::view( std::string caption ) {
    if (caption!="")
      std::cout << format("{}:\n",caption);
    auto out = this->data();
    auto m = this->m(), n = this->n();
    auto b = this->border();
    for ( idxint i=0; i<m+2*b; i++ ) {
      for ( idxint j=0; j<n+2*b; j++ ) {
        char c = ( j<n+2*b-1 ? ' ' : '\n' );
        std::cout << std::format("{:5.2}{}",out[ this->oindex(i,j) ],c);
      }
    }
  };

};

namespace sparsealg {
  template class bordered_array_1d<float>;
  template class bordered_array_1d<double>;
};
