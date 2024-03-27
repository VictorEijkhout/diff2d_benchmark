/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** linalg.cpp : bordered vector routines for omp
 ****
 ****************************************************************/

#include <cmath>
#include <algorithm>
using std::for_each;
#include <string>
using std::string;

#include <iostream>
using std::cout;
#include <format>
using std::format;

#include "../linalg.hpp"

namespace linalg {

  //! The constructor copies arguments and allocates the data
  template< typename real >
  bordered_array<real>::bordered_array( size_t m,size_t n,int border )
    : _m(m),_n(n),border(static_cast<size_t>(border))
    , std::vector<real>((m+2*border)*(n*2*border)) {}

  //! Compute the 5-point Laplace stencil from an input array
  //codesnippet d2d5ptrng
  template< typename real >
  void bordered_array<real>::central_difference_from
      ( const linalg::bordered_array<real>& other,bool trace ) {
    auto out = this->data2d();
    auto in = other.data2d();
#   pragma omp parallel for 
    for ( auto ij : inner() ) {
      auto [i,j] = ij;
      out[ i,j ] = 4*in[ i,j ]
	- in[ i-1,j ] - in[ i+1,j ] - in[ i,j-1 ] - in[ i,j+1 ];
    }
  };
  //codesnippet end

  //! Copy the interior of another bordered array, but leave border alone
  //codesnippet d2dcopyrng
  template< typename real >
  void bordered_array<real>::copy_interior_from
      ( const linalg::bordered_array<real>& other ) {
    auto out = this->data2d();
    auto in = other.data2d();
#   pragma omp parallel for 
    for ( auto ij : inner() ) {
      auto [i,j] = ij;
      out[ i,j ] = in[ i,j ];
    }
  };

  //! Scale the interior, leaving the border alone
  //codesnippet d2dscalerng
  template< typename real >
  void bordered_array<real>::scale_interior
      ( const linalg::bordered_array<real>& other, real factor ) {
    auto out = this->data2d();
    auto in = other.data2d();
#   pragma omp parallel for
    for ( auto ij : inner() ) {
      auto [i,j] = ij;
      out[ i,j ] = in[ i,j] * factor;
    }
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormrng
  template< typename real >
  real bordered_array<real>::l2norm() {
    real sum_of_squares{0};
    auto array = this->data2d();
#   pragma omp parallel for reduction(+:sum_of_squares)
    for ( auto ij : inner() ) {
      auto [i,j] = ij;
      auto v = array[i,j];
      sum_of_squares += v*v;
    }
    return std::sqrt(sum_of_squares);
  };
  //codesnippet end

  //! Set the interior to a value
  template< typename real >
  void bordered_array<real>::set( real value, bool trace ) {
    auto out = this->data2d();
    for ( auto ij : inner() ) {
      auto [i,j] = ij;
      out[i,j] = value;
    }
  };

  template< typename real >
  void bordered_array<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data2d();
#   pragma omp parallel for 
    for ( int i=0; i<_m; i++ )
      for ( int j=0; j<_n; j++ )
	if ( i==_m-1 or j==_n-1 )
	  out[ i,j ] = 1.;
  };

  template< typename real >
  void bordered_array<real>::view( string caption ) {
    if (caption!="")
      cout << format("{}:\n",caption);
    for_each
      ( this->domain().begin(),this->domain().end(),
        [n=this->_n,array = this->data2d(),border=this->border] ( auto idx ) {
          auto [i,j] = idx;
          char c = ( j<n+2*border-1 ? ' ' : '\n' );
          cout << format("{:5.2}{}",array[i,j],c);
        }
        );
  };
};

namespace linalg {
  template class bordered_array<float>;
  template class bordered_array<double>;
};
