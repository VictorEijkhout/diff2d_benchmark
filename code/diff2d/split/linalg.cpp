/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** linalg.cpp : bordered vector routines for omp
 ****    using mdspan but 1D -> 2D index conversion
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

#define HAVE_SPAN
#include "../linalg.hpp"

namespace linalg {

  //! The constructor copies arguments and allocates the data
  template< typename real >
  bordered_array<real>::bordered_array( size_t m,size_t n,int border )
    : _m(m),_n(n),border(border)
    , std::vector<real>((m+2*border)*(n*2*border)) {}

  //! Compute the 5-point Laplace stencil from an input array
  template< typename real >
  //codesnippet d2d5ptspn
  void bordered_array<real>::central_difference_from
      ( const linalg::bordered_array<real>& other,bool trace ) {
    auto out = this->data2d();
    const auto& in = other.data2d();
#   pragma omp parallel for
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->inner_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
      out[ i,j ] = 4*in[ i,j ]
	- in[ i-1,j ] - in[ i+1,j ] - in[ i,j-1 ] - in[ i,j+1 ];
    }
  };
  //codesnippet end

  //! Copy the interior of another bordered array, but leave border alone
  //codesnippet d2dcopyspn
  template< typename real >
  void bordered_array<real>::copy_interior_from
      ( const linalg::bordered_array<real>& other ) {
    auto out = this->data2d();
    auto in = other.data2d();
#   pragma omp parallel for
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->inner_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
      out[ i,j ] = in[ i,j ];
    }
  };
  //codesnippet end

  //! Scale the interior, leaving the border alone
  //codesnippet d2dscalespn
  template< typename real >
  void bordered_array<real>::scale_interior
      ( const linalg::bordered_array<real>& other, real factor ) {
    auto out = this->data2d();
    auto in = other.data2d();
#   pragma omp parallel for 
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->inner_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
      out[ i,j ] = in[ i,j ] * factor;
    }
  };
  //codesnippet end

  //! Compute the L2 norm of the interior
  //codesnippet d2dnormspn
  template< typename real >
  real bordered_array<real>::l2norm() {
    auto out = this->data2d();
    real sum_of_squares{0};
#   pragma omp parallel for reduction(+:sum_of_squares)
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->inner_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
      auto v = out[ i,j ];
      sum_of_squares += v*v;
    }
    return std::sqrt(sum_of_squares);
  };
  //codesnippet end

  //! Set the interior to a value
  template< typename real >
  void bordered_array<real>::set( real value, bool trace ) {
    auto out = this->data2d();
#   pragma omp parallel for 
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->inner_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
        out[i,j] = value;
    };
  };
  
  template< typename real >
  void bordered_array<real>::set_bc( bool down,bool right, bool trace ) {
    auto out = this->data2d();
    //#   pragma omp parallel for
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->domain_size()) ) { //this->inner() ) {
      auto [i,j] = split_i_j(ij);
      if ( i==_m-1 or j==_n-1 )
	out[ i,j ] = 1.;
    }
  };

  template< typename real >
  void bordered_array<real>::view( string caption ) {
    if (caption!="")
      cout << format("{}:\n",caption);
    auto out = this->data2d();
    for ( auto ij : rng::views::iota(static_cast<size_t>(0),this->domain_size()) ) {
      auto [i,j] = split_i_j(ij);
      char c = ( j<_n+2*border-1 ? ' ' : '\n' );
      cout << format("{:5.2}{}",out[i,j],c);
    }
  };
};

namespace linalg {
  template class bordered_array<float>;
  template class bordered_array<double>;
};
