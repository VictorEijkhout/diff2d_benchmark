/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** span.hpp : headers openmp iteration over index space
 ****
 ****************************************************************/

#ifndef SPARSEALG_SPAN_H
#define SPARSEALG_SPAN_H

#include "base.hpp"

namespace sparsealg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_span : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_span( idxint m,idxint n,int border );
#if 0
    bordered_array_span( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};
#endif
    
    void central_difference_from
        ( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;
  };

};

#endif
