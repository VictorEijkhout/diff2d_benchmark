/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** base.cpp : base routines
 ****
 ****************************************************************/

#include "base.hpp"
#include <iostream>
#include <sstream>
#include <format>

namespace sparsealg {

  /*! The constructor copies arguments and allocates the data
   * Note that the data is not zero-initialized, 
   * because this is done mode-dependent
   */
    //codesnippet d2dspan1
  template< typename real >
  bordered_array_base<real>::bordered_array_base
        ( idxint m,idxint n,int border )
    : m_(m),n_(n),_border(border)
    , data_( new real[ (m+2*border)*(n+2*border) ] )
    , data_owned(true)
    , cartesian_data
      ( md::mdspan
        ( data_,md::extents{m+2*border,n+2*border} )
        ) 
    //codesnippet end
      {
        //codesnippet d2dinner
        const auto& s = data2d();
        int b = this->border();
        idxint
          lo_m = static_cast<idxint>(b),
          hi_m = static_cast<idxint>(s.extent(0)-b),
          lo_n = static_cast<idxint>(b),
          hi_n = static_cast<idxint>(s.extent(1)-b);
        range2d = rng::views::cartesian_product
          ( rng::views::iota(lo_m,hi_m),rng::views::iota(lo_n,hi_n) );
        //codesnippet end
      };

  //! Constructor from data. This uses a zero border.
  template< typename real >
  bordered_array_base<real>::bordered_array_base( idxint m,idxint n,real *data )
    : m_(m),n_(n),_border(0)
    , data_(data)
    , data_owned(false)
    , cartesian_data( md::mdspan( data_,md::extents{m,n} ) ) 
  {};

  template< typename real >
  void bordered_array_base<real>::view( std::string caption ) {
    std::stringstream cout;
    if (caption!="")
      cout << std::format("{}:\n",caption);
    auto out = this->data();
    // auto m = this->m(), n = this->n(), n2b = this->n2b();
    // auto b = this->border();
    auto [m,n,b] = outer_sizes();
    for ( idxint i=0; i<m; i++ ) {
      for ( idxint j=0; j<n; j++ ) {
        char c = ( j<n-1 ? ' ' : '\n' );
        cout << std::format("{:5.2}{}",out[ oindex(i,j) ],c);
      }
    }
    std::cout << cout.str() << '\n';
  };

};

namespace sparsealg {
  template class bordered_array_base<float>;
  template class bordered_array_base<double>;
};
