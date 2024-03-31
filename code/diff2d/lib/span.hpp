/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** span.hpp : headers openmp iteration over index space
 ****
 ****************************************************************/

#ifndef LINALG_SPAN_H
#define LINALG_SPAN_H

#include "base.hpp"

namespace linalg {

  //codesnippet d2dbordered
  template< typename real >
  class bordered_array_span : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;
    //codesnippet end

    // constructors
    bordered_array_span( size_t m,size_t n,int border );
    bordered_array_span( size_t m,size_t n,real *data )
      : bordered_array_base<real>(m,n,data) {};

    void central_difference_from( const linalg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const linalg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set( real value,bool trace=false ) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    void view( std::string="" ) override;
  };

};

#endif
