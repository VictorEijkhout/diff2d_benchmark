#include "base.hpp"

namespace linalg {

  //! The constructor copies arguments and allocates the data
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

};

namespace linalg {
  template class bordered_array_base<float>;
  template class bordered_array_base<double>;
};
