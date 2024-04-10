/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** base.cpp : base routines
 ****
 ****************************************************************/

#include "base.hpp"
#include <iostream>
#include <sstream>
#include <format>

namespace linalg {

  /*! The constructor copies arguments and allocates the data
   * Note that the data is not zero-initialized, 
   * because this is done mode-dependent
   */
    //codesnippet d2dspan1
  template< typename real >
  bordered_array_base<real>::bordered_array_base( size_t m,size_t n,int border )
    : _m(m),_n(n),_border(static_cast<size_t>(border))
    , _n2b( n+2*border )
    , _data( new real[ (m+2*border)*(n+2*border) ] )
    , data_owned(true)
    , cartesian_data
      ( md::mdspan
        ( _data,md::extents{m+2*border,n+2*border} )
	) 
    //codesnippet end
      {};

  //! Constructor from data. This uses a zero border.
  template< typename real >
  bordered_array_base<real>::bordered_array_base( size_t m,size_t n,real *data )
    : _m(m),_n(n),_border(static_cast<size_t>(0))
    , _data(data)
    , data_owned(false)
    , cartesian_data( md::mdspan( _data,md::extents{m,n} ) ) 
  {};

  template< typename real >
  void bordered_array_base<real>::view( std::string caption ) {
    std::stringstream cout;
    if (caption!="")
      cout << std::format("{}:\n",caption);
    auto out = this->data();
    auto m = this->m(), n = this->n(), b = this->border(), n2b = this->n2b();
    for ( size_t i=0; i<m+2*b; i++ ) {
      for ( size_t j=0; j<n+2*b; j++ ) {
        char c = ( j<n+2*b-1 ? ' ' : '\n' );
        cout << std::format("{:5.2}{}",out[ oindex(i,j) ],c);
      }
    }
    std::cout << cout.str() << '\n';
  };

};

namespace linalg {
  template class bordered_array_base<float>;
  template class bordered_array_base<double>;
};
