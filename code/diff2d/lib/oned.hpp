/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** oned.hpp : headers for one-dimensional indexing for diff2d
 ****
 ****************************************************************/

#ifndef SPARSEALG_ONED_H
#define SPARSEALG_ONED_H

#include <tuple>
#include <vector>

#include "base.hpp"

namespace sparsealg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_1d : public bordered_array_base<real> {
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_1d( idxint m,idxint n,int border );
    bordered_array_1d( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};

    // required functionality
    void central_difference_from
        ( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;

    /*
     * utility section
     */
    
    //! convert linear ij interior to pair i,j in the global data
    std::pair<idxint,idxint> split_i_j( idxint ij ) {
      return std::make_pair( ij/this->n()+this->border(), ij%this->n()+this->border() );
    };

    //! m/n/border size of the allocated data
    std::tuple<idxint,idxint,idxint> outer_sizes() { 
      return std::make_tuple
        ( this->m()+2*this->border(),
          this->n()+2*this->border(),this->border() ); };

    //! m/n size of the domain, and the border
    std::tuple<idxint,idxint,idxint> inner_sizes() { 
      return std::make_tuple( this->m(),this->n(),this->border() ); };

    std::vector<real> internal_data();

  };

};

#endif
