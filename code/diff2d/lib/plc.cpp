/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** plc.cpp : range policy based implementation
 ****
 ****************************************************************/

#include <cmath>
#include <algorithm>
#include <execution>
#include <string>

#include <iostream>
using std::cout;
#include <format>
using std::format;

#include "range.hpp"

namespace linalg {

  template <typename... _IteratorTypes>
  using __are_random_access_iterators
  = std::__and_<std::is_base_of<std::random_access_iterator_tag,
				std::__iterator_category_t<_IteratorTypes>>...>;

    // range algorithms on ranges doesn't work
    // auto inner = this->inner();
  template< typename real >
  void bordered_array_range<real>::central_difference_from
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {
    const auto& other = dynamic_cast<const linalg::bordered_array_range<real>&>(_other);
    auto& out = this->data2d();
    const auto& in = other.data2d();
    auto inner = inner_range<real>(*this);
    std::for_each
      ( std::execution::par,
        inner.begin(),inner.end(),
        [out,in] ( auto idx ) {
          auto [i,j] = idx;
          out[i,j] = 4*in[i,j] - in[i-1,j] - in[i+1,j] - in[i,j-1] - in[i,j+1];
        }
      );
  };

  template< typename real >
  void bordered_array_range<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
    const auto& other = dynamic_cast<const linalg::bordered_array_range<real>&>(_other);
    auto& out = this->data2d();
    const auto& in = other.data2d();
    auto inner = inner_range<real>(*this);
    std::for_each
      ( std::execution::par,
        inner.begin(),inner.end(),
        [out,in,factor] ( auto idx ) {
          auto [i,j] = idx;
          out[i,j] = in[i,j] * factor;
      }
    );
  };

  template< typename real >
  real bordered_array_range<real>::l2norm() {
    real sum_of_squares{0};
    const auto& array = this->data2d();
    auto inner = this->inner();
    std::for_each
      ( std::execution::par,
        inner.begin(),inner.end(),
        [array,&sum_of_squares] ( auto idx ) {
          auto [i,j] = idx;
          sum_of_squares += array[i,j] * array[i,j];
      }
    );
    return std::sqrt(sum_of_squares);
  };

  //! Set the interior to a value
  template< typename real >
  void bordered_array_range<real>::set( real value, bool trace ) {
    auto& out = this->data2d();
    std::for_each
      ( std::execution::par,
        this->domain().begin(),this->domain().end(),
        [out,value] ( auto idx ) {
          auto [i,j] = idx;
          out[i,j] = value;
      }
    );
  };

  template< typename real >
  void bordered_array_range<real>::set_bc( bool down,bool right, bool trace ) {
    auto& out = this->data2d();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    std::for_each
      ( std::execution::par,
        this->domain().begin(),this->domain().end(),
        [m,n,out] ( auto idx ) {
          auto [i,j] = idx;
          if ( i==m-1 or j==n-1 )
            out[ i,j ] = 1.;
        }
      );
  };

  template< typename real >
  void bordered_array_range<real>::view( std::string caption ) {
    if (caption!="")
      std::cout << format("{}:\n",caption);
    auto& out = this->data2d();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    auto inner = this->inner();
    std::for_each
      ( std::execution::par,
        inner.begin(),inner.end(),
        [m,n,n2b,out] ( auto idx ) {
        auto [i,j] = idx;
        if ( i==m-1 or j==n-1 ) {
          char c = ( j<n2b-1 ? ' ' : '\n' );
          cout << format("{:5.2}{}",out[i,j],c);
        } }
        );
  };

};

namespace linalg {
  template class bordered_array_range<float>;
  template class bordered_array_range<double>;
};
