/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diy2d.hpp : based on span.hpp but my own cartesian rangee
 ****
 ****************************************************************/

#ifndef SPARSEALG_DIY2D_H
#define SPARSEALG_DIY2D_H

#include "base.hpp"

namespace sparsealg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_diy2d : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_diy2d( idxint m,idxint n,int border );
#if 0
    bordered_array_diy2d( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};
#endif
    
    void central_difference_from( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;

    /*
     * Iterator class
     */
    class cartesian_product_diy {
    private: 
      idxint m, n; int b;
    public:
      class cartesian_iterator {
      private:
        idxint m, n; int b; /* global domain description */
        idxint i{0}, j{0};  /* local iteration point */
      public:
        cartesian_iterator\
            ( idxint m, idxint n, int b, idxint i, idxint j )
              : m(m),n(n),b(b),i(i),j(j) {};
        bool operator==( const cartesian_iterator& other ) const {
          return this->j==other.j and this->i==other.i; };
        auto operator*() const {
          return std::make_pair( (i+b), (j+b) ); };
        //codesnippet d2ddiyiter
        auto& operator++(  ) {
          j++; i+= (j/m); j = j%m;
          return *this; };
        //codesnippet end
        auto operator+( idxint dist ) const {
          auto displaced(*this);
          auto lin = ( i*m+j ) + dist;
          displaced.i = lin/m; displaced.j = lin%m; 
          return displaced;
        };
        // gcc needs +=, intel only +
        auto& operator+=( idxint dist ) {
          auto lin = ( i*m+j ) + dist;
          i = lin/m; j = lin%m; 
          return *this;
        };
        idxint operator-( const cartesian_iterator& other ) const {
          return (i-other.i)*m + (j-other.j); }
      };
      cartesian_product_diy( idxint m, idxint n,int b ) 
        : m(m),n(n),b(b) {};
      auto begin() {
        return cartesian_iterator(m,n,b, /* start */ 0,0); };
      auto end()   {
        return cartesian_iterator(m,n,b, /* end   */ m,0); };
    };
    auto inner_diy() { 
      return cartesian_product_diy(this->m(),this->n(),this->border());
    };
  };
};

#endif
