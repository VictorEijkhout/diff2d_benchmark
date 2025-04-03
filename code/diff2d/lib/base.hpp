/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** base.hpp : headers for pure virtual base classes for diff2d implementation
 ****
 ****************************************************************/

#ifndef SPARSEALG_BASE_H
#define SPARSEALG_BASE_H

#include <cstddef>
#include <cstdint>
using idxint = std::int64_t;
using uidxint = std::uint64_t;

#include <memory>
#include <optional>
#include <string>
#include <tuple>

#include <ranges>
namespace rng = std::ranges;
#include <tuple>
#include <vector>

#include "mdspan/mdspan.hpp"
namespace md = Kokkos;

class logging {
private:
  static inline float flopcount{0},bytecount{0};
public:
  void log_flops( float n ) { flopcount += n; };
  void log_bytes( float n ) { bytecount += n; };
  std::pair<float,float> log_report() {
    return std::make_pair(flopcount,bytecount); };
};

//codesnippet d2dseqindexf
inline idxint index2d( idxint i,idxint j,idxint m,idxint n,idxint b ) {
  return (i+b)*(n+2*b) + j+b;
};
//codesnippet end

namespace sparsealg {

  // forward definition so that we can declare it `friend'
  template< typename real >
  class distributed_array;

  //codesnippet d2dvirtualbase
  template< typename real >
  class bordered_array_base {
    //codesnippet end
    friend class distributed_array<real>;
  protected:
    idxint m_{0},n_{0}; int _border{0}; 
  public:
    virtual ~bordered_array_base() {
      // if (data_owned) { delete[] data_; }
    };
    //codesnippet d2dbaseconstruct
    bordered_array_base( idxint m,idxint n,int border );
    //codesnippet end
    //! We need the default constructor for the MPI version
    bordered_array_base() = default;
    // Constructor from data, only for the MPI gathered version
    bordered_array_base( idxint m,idxint n,real *data );
    std::tuple<idxint,idxint,int> outer_sizes() const {
      return std::make_tuple(m_+2*_border,n_+2*_border,_border); };
    std::tuple<idxint,idxint,int,idxint,idxint> inner_sizes() const {
      return std::make_tuple(m_,n_,_border,m_+2*_border,n_+2*_border); };
    inline auto m() const { return m_; };
    inline auto n() const { return n_; };
    inline auto border() const { return _border; };
    inline auto n2b() const { return n_+2*_border; };
  private:
    logging _logger;

    //codesnippet d2dspan0
  private:
    std::unique_ptr<real[]> data_{nullptr};
    md::mdspan<
      real,
      md::dextents<idxint,2>
              > cartesian_data;
    //codesnippet end
    bool data_owned;

  protected:
    /*
     * Data access
     */

    //codesnippet d2dbaredata
    real* data() { return data_.get(); };
    real const * const data() const { return data_.get(); };
    //codesnippet end

    //codesnippet d2dspan2
    //! pointer to the data as 2D array
    auto& data2d() {
      return cartesian_data; };
    const auto& data2d() const {
      return cartesian_data; };
    //codesnippet end

    // template<typename Self>
    // auto& data2d( this Self& self ) {
    //   return self.cartesian_data;
    // };

    //! index in the bordered array
    auto oindex( int i,int j ) {
      return i*(this->n()+2*this->border()) + j; };
    //! index in the interior
    auto iindex( int i,int j ) {
      return (i+this->border())*(this->n()+2*this->border()) + j+this->border(); };

    /*! A view to the interior
     * \todo we should really use m_,n_ instead of the extents
     * \todo we really want cbegin / cend but not yet available
     */
     mutable std::optional< decltype( rng::views::cartesian_product
                             ( rng::views::iota(idxint{0},idxint{0}),
                               rng::views::iota(idxint{0},idxint{0}) ) ) >
            range2d_ = {};
  private:
    decltype( rng::views::cartesian_product
	      ( rng::views::iota(idxint{0},idxint{0}),
		rng::views::iota(idxint{0},idxint{0}) ) )  range2d;
  public:
    //codesnippet d2dinnerrange
    auto& inner_range() const { return range2d; };
    //codesnippet end
#if 0
      if (not range2d.has_value()) {
        //codesnipper d2dinner
        const auto& s = data2d();
        int b = this->border();
        idxint
          lo_m = static_cast<idxint>(b),
          hi_m = static_cast<idxint>(s.extent(0)-b),
          lo_n = static_cast<idxint>(b),
          hi_n = static_cast<idxint>(s.extent(1)-b);
        range2d = rng::views::cartesian_product
          ( rng::views::iota(lo_m,hi_m),rng::views::iota(lo_n,hi_n) );
        //codesnipper end
      }
      return *range2d;
    };
#endif

    mutable std::optional< decltype( rng::views::iota(idxint{0},idxint{0}) ) >
            range2di = {};
    auto inneri() const {
      if (not range2di.has_value()) {
        const auto& s = data2d();
        int b = this->border();
        idxint
          lo_m = static_cast<idxint>(b),
          hi_m = static_cast<idxint>(s.extent(0)-b),
          lo_n = static_cast<idxint>(b),
          hi_n = static_cast<idxint>(s.extent(1)-b);
        range2di = rng::views::iota(lo_m,hi_m);
      }
      return *range2di;
    };

    mutable std::optional< decltype( rng::views::iota(idxint{0},idxint{0}) ) >
            range2dj = {};
    auto innerj() const {
      if (not range2dj.has_value()) {
        const auto& s = data2d();
        int b = this->border();
        idxint
          lo_m = static_cast<idxint>(b),
          hi_m = static_cast<idxint>(s.extent(0)-b),
          lo_n = static_cast<idxint>(b),
          hi_n = static_cast<idxint>(s.extent(1)-b);
        range2dj = rng::views::iota(lo_n,hi_n);
      }
      return *range2dj;
    };

    auto inner_size() const {
      const auto& s = data2d();
      int b = this->border();
      idxint
        lo_m = static_cast<idxint>(b),
        hi_m = static_cast<idxint>(s.extent(0)-b),
        lo_n = static_cast<idxint>(b),
        hi_n = static_cast<idxint>(s.extent(1)-b);
      return (hi_m-lo_m) * ( hi_n-lo_n);
    };

    auto domain() {
      const auto& s = data2d();
      return rng::views::cartesian_product
        ( rng::views::iota(0,static_cast<int>( s.extent(0) )),
          rng::views::iota(0,static_cast<int>( s.extent(1) ))
          );
    };

  public:
    // essential functionality
    //codesnippet d2dvirtualfunc
    virtual void central_difference_from
      ( const sparsealg::bordered_array_base<real>&,bool=false ) = 0;
    virtual void scale_interior
      ( const sparsealg::bordered_array_base<real>&, real ) = 0;
    virtual real l2norm() = 0;
    virtual void set_value( real value,bool trace=false ) = 0;
    virtual void set_bc(bool down, bool right, bool trace=false) = 0;
    //codesnippet end
    virtual void view( std::string="" );

    // logging
    void log_flops( float n ) { _logger.log_flops(n); };
    void log_bytes( float n ) { _logger.log_bytes(n); };
    std::pair<float,float> log_report() { return _logger.log_report(); };
  };

  template< typename real >
  class inner_range {
  private:
    const bordered_array_base<real> &a;
  public:
    inner_range( const bordered_array_base<real>& a ) : a(a) {};
    class iter {
    private:
      idxint lo_m,hi_m, lo_n,hi_n;
      int border;
      idxint n2b;
      idxint seek_i,seek_j;
    public:
      using coordinate = std::pair<idxint,idxint>;
      using iterator_category = std::random_access_iterator_tag;
      using value_type = coordinate;
      using difference_type = idxint; // std::ptrdiff_t; // integral type
      using iter_difference_type = idxint; // std::ptrdiff_t; // integral type

      using pointer = coordinate*;
      using pointer_type = coordinate*;
      using reference  = coordinate&;
      using reference_type  = coordinate&;
      // constructor
      iter() = default;
      iter( idxint lo_m,idxint hi_m, idxint lo_n,idxint hi_n,
            int border, idxint n2b, idxint seek_i,idxint seek_j )
        : lo_m(lo_m),hi_m(hi_m), lo_n(lo_n),hi_n(hi_n)
        , border(border), n2b(n2b)
        , seek_i(seek_i), seek_j(seek_j) {};
      // operators
      auto operator*() const { return std::make_pair(seek_i,seek_j); };

      auto operator==( const iter& other ) const {
        return seek_i==other.seek_i or seek_j==other.seek_j; };
      auto operator!=( const iter& other ) const {
        return seek_i!=other.seek_i or seek_j!=other.seek_j; };

      auto operator<( const iter& other ) const {
        return seek_i<other.seek_i or
                      ( seek_i==other.seek_i and seek_j<other.seek_j ); }
      auto operator>( const iter& other ) const {
        return seek_i>other.seek_i or
                      ( seek_i==other.seek_i and seek_j>other.seek_j ); }
      auto operator>=( const iter& other ) const {
        return operator>(other) or operator==(other); };
      auto operator<=( const iter& other ) const {
        return operator<(other) or operator==(other); };

      auto& operator++() { // needs to return by reference!
        seek_j++; if (seek_j==hi_n) { seek_j=lo_n; seek_i++; }
        return *this; };
      auto& operator--() { // needs to return by reference!
        seek_j--; if (seek_j==border) { seek_j=hi_n; seek_i--; }
        return *this; };
      auto operator++(int) { auto tmp(*this); ++(*this); return tmp; };
      auto operator--(int) { auto tmp(*this); --(*this); return tmp; };

      idxint operator-( const iter& other ) const {
        auto lin1 = (seek_i-border) * ( hi_n-lo_n ) + seek_j-border;
        auto lin2 = (other.seek_i-other.border) * ( hi_n-lo_n ) + other.seek_j-other.border;
        return lin1-lin2; };
      auto& operator+=( idxint dist ) {
        auto lines = dist / (hi_n-lo_n);
        auto pts = dist - lines*(hi_n-lo_n);
        seek_i += lines; seek_j += pts;
        if ( seek_j> hi_n ) { seek_i++; seek_j -= (hi_n-lo_n); };
        return *this; };
      auto operator+( idxint dist ) const {
        auto tmp(*this); tmp += dist; return tmp; };

    };

    auto begin() { return iter
        ( a.border(),a.border()+a.m(), a.border(),a.border()+a.n(),
          a.n2b(), a.border(),
          a.border(),a.border() ); };
    auto end() { return iter
        ( a.border(),a.border()+a.m(), a.border(),a.border()+a.n(),
          a.n2b(), a.border(),
          a.border()+a.m(),a.border() ); };
  };

  template< typename real >
  auto operator+( idxint dist,const typename inner_range<real>::iter& cur ) {
    return cur+dist; };

  // static_assert( std::input_iterator< inner_range<float>::iter > );
  // static_assert( std::sentinel_for< inner_range<float>::iter,inner_range<float>::iter > );
  // static_assert( std::input_iterator< inner_range<double>::iter > );
  // static_assert( std::sentinel_for< inner_range<double>::iter,inner_range<double>::iter > );

  // static_assert( std::ranges::random_access_range<inner_range<float>> );
  // static_assert( std::ranges::random_access_range<inner_range<double>> );
};

#endif

