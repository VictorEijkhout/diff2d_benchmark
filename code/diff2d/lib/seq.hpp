/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** seq.hpp : headers for one-dimensional indexing for diff2d
 ****
 ****************************************************************/

#ifndef SPARSEALG_SEQ_H
#define SPARSEALG_SEQ_H

#include <tuple>
#include <vector>

#include "base.hpp"

namespace sparsealg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_seq : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_seq( idxint m,idxint n,int border );
#if 0
    bordered_array_seq( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};
#endif

    // required functionality
    void central_difference_from( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;

    /*
     * utility section
     */
    
    //! convert linear ij interior to pair i,j in the global data
    std::pair<idxint,idxint> split_i_j( idxint ij ) {
      return std::make_pair( ij/this->n()+this->border(), ij%this->n()+this->border() );
    };

    //! index in the bordered array
    inline auto oindex( int i,int j ) {
      return i*(this->n()+2*this->border()) + j; };
    //! index in the interior
    inline auto iindex( int i,int j ) {
      return (i+this->border())*(this->n()+2*this->border()) + j+this->border(); };

    std::vector<real> internal_data();

  };

};

#endif
