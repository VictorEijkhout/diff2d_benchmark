/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** iota.hpp : headers openmp iteration over index space
 ****
 ****************************************************************/

#ifndef SPARSEALG_IOTA_H
#define SPARSEALG_IOTA_H

#include "base.hpp"

namespace sparsealg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_iota : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_iota( idxint m,idxint n,int border );
    bordered_array_iota( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};

    void central_difference_from( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;
  };

};

#endif
