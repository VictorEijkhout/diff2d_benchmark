/****************************************************************
 *  ****
 *  **** This file belongs with the course
 *  **** Parallel Programming in MPI and OpenMP
 *  **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 *  ****
 *  **** sycl.hpp : headers for sycl implementation for diff2d
 *  **** code contributed by Yojan Chitkara
 *  ****
 *  ****************************************************************/

#ifndef SPARSEALG_SYCL_H
#define SPARSEALG_SYCL_H

#include <tuple>
#include <vector>

#include <sycl/sycl.hpp>
using namespace sycl;

#include "base.hpp"


namespace sparsealg {

  //codesnippet d2dborderedsycl
  template< typename real >
  class bordered_array_sycl : public bordered_array_base<real> {
    friend class distributed_array<real>;
  public:
    sycl::queue q;
    bordered_array_sycl( idxint m,idxint n,int border,sycl::queue q );
    //codesnippet end
    bordered_array_sycl( idxint m,idxint n,real *data )
      : bordered_array_base<real>(m,n,data) {};

    // bookkeeping
    using bordered_array_base<real>::log_flops;
    using bordered_array_base<real>::log_bytes;

    //required functionality
    void central_difference_from( const sparsealg::bordered_array_base<real>&,bool=false ) override;
    void scale_interior( const sparsealg::bordered_array_base<real>&, real ) override;
    real l2norm() override;
    void set_value( real value,bool trace=false ) override;
    void view(std::string caption) override;
    void set_bc(bool down, bool right, bool trace=false) override;
    std::vector<real> internal_data();
  };
};

#endif
