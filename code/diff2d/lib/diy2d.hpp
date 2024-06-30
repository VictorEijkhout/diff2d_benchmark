/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diy2d.hpp : based on span.hpp but my own cartesian rangee
 ****
 ****************************************************************/

#ifndef LINALG_DIY2D_H
#define LINALG_DIY2D_H

#include "base.hpp"

namespace linalg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_diy2d : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_diy2d( int64_t m,int64_t n,int border );
    bordered_array_diy2d( int64_t m,int64_t n,real *data )
      : bordered_array_base<real>(m,n,data) {};

    void central_difference_from( const linalg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const linalg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;
  };

  /*
   * Iterator class
   */
  class twod_iterator {
  private: 
    std::int64_t m, n; int b;
    std::int64_t i{0}, j{0}; 
  public:
    two_iterator( std::int64_t m, std::int64_t n,int b ) 
      : m(m),n(n),b(b) {};
    auto& begin() { return *this; };
    auto& end() { i=m; j=0; return *this; };
    bool operator!=( const two_iterator& other ) const {
      return this->j!=other->j or this->i!=other->i; };
    std::int64_t operator*() const {
      return (i+b)*(m+2*b) + (j+b); };
    auto& operator++( int ) {
      j++; i+= (j/m); j = j%m; return *this; };
  };
  auto inner() override {
    return twod_iterator(this->m(),this->n(),this->border());
  };

};

#endif
