/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** base.hpp : headers for pure virtual base classes for diff2d implementation
 ****
 ****************************************************************/

#ifndef LINALG_BASE_H
#define LINALG_BASE_H

#include <cstddef>
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

namespace linalg {

  // forward definition so that we can declare it `friend'
  template< typename real >
  class distributed_array;

  //codesnippet d2dvirtualbase
  template< typename real >
  class bordered_array_base {
    //codesnippet end
    friend class distributed_array<real>;
  protected:
    std::int64_t _m{0},_n{0}; int _border{0}; 
  public:
    virtual ~bordered_array_base() {
      if (data_owned) { delete[] _data; } };
    //codesnippet d2dbaseconstruct
    bordered_array_base( std::int64_t m,std::int64_t n,int border );
    //codesnippet end
    //! We need the default constructor for the MPI version
    bordered_array_base() = default;
    // Constructor from data, only for the MPI gathered version
    bordered_array_base( std::int64_t m,std::int64_t n,real *data );
    std::tuple<std::int64_t,std::int64_t,int> outer_sizes() const {
      return std::make_tuple(_m+2*_border,_n+2*_border,_border); };
    std::tuple<std::int64_t,std::int64_t,int,std::int64_t,std::int64_t> inner_sizes() const {
      return std::make_tuple(_m,_n,_border,_m+2*_border,_n+2*_border); };
    inline auto m() const { return _m; };
    inline auto n() const { return _n; };
    inline auto border() const { return _border; };
    inline auto n2b() const { return _n+2*_border; };
  private:
    logging _logger;

    //codesnippet d2dspan0
  private:
    real *_data{nullptr};
    bool data_owned;
    md::mdspan< real,md::dextents<std::int64_t,2> > 
        cartesian_data;
    //codesnippet end

  protected:
    /*
     * Data access
     */

    //codesnippet d2dbaredata
    real* data() { return _data; };
    real const * const data() const { return _data; };
    //codesnippet end

    //codesnippet d2dspan2
    //! pointer to the data as 2D array
    auto& data2d() {
      return cartesian_data;
    };
    //codesnippet end
    const auto& data2d() const {
      return cartesian_data;
    };

    //! index in the bordered array
    auto oindex( int i,int j ) {
      return i*(this->n()+2*this->border()) + j; };
    //! index in the interior
    auto iindex( int i,int j ) {
      return (i+this->border())*(this->n()+2*this->border()) + j+this->border(); };

    /*! A view to the interior
     * \todo we should really use _m,_n instead of the extents
     * \todo we really want cbegin / cend but not yet available
     */
    //codesnippet d2dinner
    auto inner() {
      const auto& s = data2d();
      int b = this->border();
      std::int64_t
        lo_m = static_cast<std::int64_t>(b),
        hi_m = static_cast<std::int64_t>(s.extent(0)-b),
        lo_n = static_cast<std::int64_t>(b),
        hi_n = static_cast<std::int64_t>(s.extent(1)-b);
      return rng::views::cartesian_product
        ( rng::views::iota(lo_m,hi_m),rng::views::iota(lo_n,hi_n) );
    };
    //codesnippet end

    auto inner_size() const {
      const auto& s = data2d();
      int b = this->border();
      std::int64_t
        lo_m = static_cast<std::int64_t>(b),
        hi_m = static_cast<std::int64_t>(s.extent(0)-b),
        lo_n = static_cast<std::int64_t>(b),
        hi_n = static_cast<std::int64_t>(s.extent(1)-b);
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
      ( const linalg::bordered_array_base<real>&,bool=false ) = 0;
    virtual void scale_interior( const linalg::bordered_array_base<real>&, real ) = 0;
    virtual real l2norm() = 0;
    virtual void set( real value,bool trace=false ) = 0;
    virtual void set_bc(bool down, bool right, bool trace=false) = 0;
    virtual void view( std::string="" );
    //codesnippet end

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
      std::int64_t lo_m,hi_m, lo_n,hi_n;
      int border;
      std::int64_t n2b;
      std::int64_t seek_i,seek_j;
    public:
      using coordinate = std::pair<std::int64_t,std::int64_t>;
      using iterator_category = std::random_access_iterator_tag;
      using value_type = coordinate;
      using difference_type = std::int64_t; // std::ptrdiff_t; // integral type
      using iter_difference_type = std::int64_t; // std::ptrdiff_t; // integral type

      using pointer = coordinate*;
      using pointer_type = coordinate*;
      using reference  = coordinate&;
      using reference_type  = coordinate&;
      // constructor
      iter() = default;
      iter( std::int64_t lo_m,std::int64_t hi_m, std::int64_t lo_n,std::int64_t hi_n,
            int border, std::int64_t n2b, std::int64_t seek_i,std::int64_t seek_j )
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

      std::int64_t operator-( const iter& other ) const {
        auto lin1 = (seek_i-border) * ( hi_n-lo_n ) + seek_j-border;
        auto lin2 = (other.seek_i-other.border) * ( hi_n-lo_n ) + other.seek_j-other.border;
        return lin1-lin2; };
      auto& operator+=( std::int64_t dist ) {
        auto lines = dist / (hi_n-lo_n);
        auto pts = dist - lines*(hi_n-lo_n);
        seek_i += lines; seek_j += pts;
        if ( seek_j> hi_n ) { seek_i++; seek_j -= (hi_n-lo_n); };
        return *this; };
      auto operator+( std::int64_t dist ) const {
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
  auto operator+( std::int64_t dist,const typename inner_range<real>::iter& cur ) {
    return cur+dist; };

  // static_assert( std::input_iterator< inner_range<float>::iter > );
  // static_assert( std::sentinel_for< inner_range<float>::iter,inner_range<float>::iter > );
  // static_assert( std::input_iterator< inner_range<double>::iter > );
  // static_assert( std::sentinel_for< inner_range<double>::iter,inner_range<double>::iter > );

  // static_assert( std::ranges::random_access_range<inner_range<float>> );
  // static_assert( std::ranges::random_access_range<inner_range<double>> );
};

#endif

