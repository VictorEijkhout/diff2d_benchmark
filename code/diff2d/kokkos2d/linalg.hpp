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

class logging {
private:
  static inline float flopcount{0},bytecount{0};
public:
  void log_flops( float n ) { flopcount += n; };
  void log_bytes( float n ) { bytecount += n; };
  std::pair<float,float> log_report() {
    return std::make_pair(flopcount,bytecount); };
};

namespace linalg {

  template< typename real >
  class distributed_array;

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array : public std::vector<real> {
    friend class distributed_array<real>;
  private:
    size_t _m{0},_n{0},border{0}; 
    size_t _n2b{0}; // nsize with borders: only used to optimize linear indexing
    logging _logger;
  public:
    //! We need the default constructor for the MPI version
    bordered_array() = default;
    bordered_array( size_t m,size_t n,int border=1 );
    //codesnippet end
    //! Constructor from data. This uses a zero border.
    bordered_array( size_t m,size_t n,real *data )
      : _m(m),_n(n),border(0),
	std::vector<real>( data,data+m*n ) { _n2b = _n+2*border; };

    // logging
    void log_flops( float n ) { _logger.log_flops(n); };
    void log_bytes( float n ) { _logger.log_bytes(n); };
    std::pair<float,float> log_report() { return _logger.log_report(); };

    //codesnippet d2dspan
    //! pointer to the data as 2D array
    auto data2d() {
      return md::mdspan
        ( std::vector<real>::data(),
          md::extents{_m+2*border,_n+2*border} );
    };
    //codesnippet end
    //! const pointer to the data as 2D array
    auto data2d() const {
      return md::mdspan
        ( std::vector<real>::data(),
          md::extents{_m+2*border,_n+2*border} );
    };

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
    //! convert linear ij interior to pair i,j in the global data
    std::pair<size_t,size_t> split_i_j( size_t ij ) {
      return std::make_pair( ij/_n+border, ij%_n+border );
    };

    //! index in the bordered array
    inline auto oindex( int i,int j ) { return i*(_n+2*border) + j; };
    //! index in the interior
    inline auto iindex( int i,int j ) { return (i+border)*(_n+2*border) + j+border; };

    //! m/n/border size of the allocated data
    std::tuple<size_t,size_t,size_t> outer_sizes() { 
      return std::make_tuple( _m+2*border,_n+2*border,border ); };

    //! m/n size of the domain, and the border
    std::tuple<size_t,size_t,size_t> inner_sizes() { 
      return std::make_tuple( _m,_n,border ); };

    std::vector<real> internal_data();

    void central_difference_from( const linalg::bordered_array<real>&,bool=false );
    void copy_interior_from( const linalg::bordered_array<real>& );
    void scale_interior( const linalg::bordered_array<real>&, real );
    real l2norm();
    void set( real value,bool trace=false );
    void set_bc(bool down, bool right, bool trace=false);
    void view( std::string="" );
  };
};

#endif
