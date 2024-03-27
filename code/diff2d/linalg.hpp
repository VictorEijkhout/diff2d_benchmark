/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2023 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** linalg.hpp : headers for bordered vector
 ****
 ****************************************************************/

#ifndef LINALG_HPP
#define LINALG_HPP

#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <ranges>
namespace rng = std::ranges;

#include "mdspan/mdspan.hpp"
namespace md = Kokkos;
namespace KokkosEx = MDSPAN_IMPL_STANDARD_NAMESPACE::MDSPAN_IMPL_PROPOSED_NAMESPACE;
namespace mdx = KokkosEx;

namespace linalg {
  
  template< typename real >
  class bordered_array_2d : public bordered_array_base<real> {
  public:
    //codesnippet d2dspan
    //! pointer to the data as 2D array
    auto data2d() {
      return md::mdspan
        ( std::vector<real>::data(),
          md::extents{_m+2*border,_n+2*border} );
    };
    //! const pointer to the data as 2D array
    auto data2d() const {
      return md::mdspan
        ( std::vector<real>::data(),
          md::extents{_m+2*border,_n+2*border} );
    };
    //codesnippet end

    //! A view the full data
    auto domain() {
      const auto& s = data2d();
      auto
	z = static_cast<size_t>( 0 ),
	m = static_cast<size_t>( s.extent(0) ),
	n = static_cast<size_t>( s.extent(1) );
      return rng::views::cartesian_product
	( rng::views::iota(z,m),rng::views::iota(z,n) );
    }

    auto domain_size() {
      return (_m+2*border)*(_n+2*border);
    };
    
    //codesnippet d2dinner
    /*! A view to the interior
     * \todo we should really use _m,_n instead of the extents
     */
    auto inner() {
      const auto s = data2d();
      // we really want cbegin / cend but not yet available
      size_t
	lo_m = static_cast<size_t>(border),
	hi_m = static_cast<size_t>(s.extent(0)-border),
	lo_n = static_cast<size_t>(border),
	hi_n = static_cast<size_t>(s.extent(1)-border);
      return rng::views::cartesian_product
	( rng::views::iota(lo_m,hi_m),rng::views::iota(lo_n,hi_n) );
    };
    //codesnippet end
    // and a const version
    const auto inner() const {
      const auto s = data2d();
      // we really want cbegin / cend but not yet available
      size_t
	lo_m = static_cast<size_t>(border),
	hi_m = static_cast<size_t>(s.extent(0)-border),
	lo_n = static_cast<size_t>(border),
	hi_n = static_cast<size_t>(s.extent(1)-border);
      return rng::views::cartesian_product
	( rng::views::iota(lo_m,hi_m),rng::views::iota(lo_n,hi_n) );
    };

    size_t inner_size() {
      const auto s = data2d();
      auto 
	sm = static_cast<size_t>( s.extent(0) )-2*border,
	sn = static_cast<size_t>( s.extent(1) )-2*border;
      return sm*sn;
    };
  };

};

#endif
