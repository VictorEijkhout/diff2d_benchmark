/****************************************************************
 *  ****
 *  **** This file belongs with the course
 *  **** Parallel Programming in MPI and OpenMP
 *  **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 *  ****
 *  **** sycl.cpp : 2D diffusion in parallel through SYCL 
 *  **** code contributed by Yojan Chitkara
 *  ****
 *****************************************************************/

#include <cmath>
#include <format>
#include <iostream>
#include <string>

#include "sycl.hpp"

namespace sparsealg {
  template< typename real >
  bordered_array_sycl<real>::bordered_array_sycl( idxint m,idxint n,int border, queue q )
    : bordered_array_base<real>(m,n,border) {
    auto out = this->data();

    q = this->q;
    buffer<real,2> Buf_a(out, sycl::range<2>(m,n));

    q.submit([&](sycl::handler &h) {
      sycl::accessor D_a(Buf_a, h, sycl::write_only);

      h.parallel_for(sycl::range<2>(m, n), [=](auto index) {
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;
        D_a[row][col] = static_cast<real>(0);
      });
    }).wait();  
  };

  template< typename real >
  void bordered_array_sycl<real>::central_difference_from
      ( const sparsealg::bordered_array_base<real>& _other,bool trace ) {

    const auto& other = dynamic_cast<const sparsealg::bordered_array_sycl<real>&>(_other);

    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    auto q = this->q;
    //codesnippet d2d5ptsycl
    sycl::range<2> outer_range(m2b,n2b);
    buffer<real,2> Buf_a(out, outer_range);
    buffer<real,2> Buf_b(in,  outer_range);
    sycl::range<2> inner_range(m,n);
    sycl::id<2> offset{1,1};
    q.submit([&] (handler &h)
    {
      accessor D_a(Buf_a,h,inner_range,offset,read_only);
      accessor D_b(Buf_b,h,inner_range,offset,write_only);

      h.parallel_for(inner_range, [=](item<2> index){
        auto row = index.get_id(0);
        auto col = index.get_id(1);

        real stencil_value =
          4*D_a[row][col]
          - D_a[row-1][col] - D_a[row+1][col]
          - D_a[row][col-1] - D_a[row][col+1];
        D_b[row][col] = stencil_value;
      });
    }).wait();
    //codesnippet end
  };

  template < typename real >
  void bordered_array_sycl<real>::set_value( real value, bool trace ) {
    auto out = this->data();
    auto [m_,n_,b,m2b_,n2b_] = this->inner_sizes();
    auto m = static_cast<uidxint>(m_);
    auto n = static_cast<uidxint>(n_);
    auto m2b = static_cast<uidxint>(m2b_);
    auto n2b = static_cast<uidxint>(n2b_);

    auto q = this->q;
    //codesnippet syclbufaccess
    buffer<real,2> Buf_a(out, sycl::range<2>(m,n));
    q.submit([&](sycl::handler &h) {
      sycl::accessor D_a(Buf_a, h, sycl::write_only);
    //codesnippet end

      h.parallel_for(sycl::range<2>(m, n),
        [=](auto index) {
          auto row = index.get_id(0)+1;
          auto col = index.get_id(1)+1;
          D_a[row][col] = value;
          });
    }).wait();

    return;
  };

  //codesnippet d2dnormsycl
  template < typename real >
  real bordered_array_sycl<real>::l2norm() {
    real norm = 0.0;
    auto out = this->data();
    auto m = std::get<0>(this->inner_sizes());
    auto n = std::get<0>(this->inner_sizes());
    buffer<real,2> Buf_a
      (out, sycl::range<2>(m,n));
    buffer<real,1> Buf_n
      (&norm, range<1>(1));

    auto q = this->q;
    q.submit([&](handler &h) {
      accessor D_a(Buf_a, h, read_only);
      auto D_Fn = reduction
	(Buf_n, h, std::plus<real>());

      h.parallel_for
        (range<2>(m - 2, n - 2),
         D_Fn,
         [=](item<2> index, auto &sum) {
           auto row = index.get_id(0) + 1;
           auto col = index.get_id(1) + 1;
           real value = D_a[row-1][col-1];
           sum += (value*value);
      });
    }).wait();
    
    norm = std::sqrt(norm);
    return norm;
  };
  //codesnippet end

  template < typename real >
  void bordered_array_sycl<real>::view( std::string caption ) {
    if (caption!="")
      std::cout << format("{}:\n",caption);

    auto out = this->data();
    auto m = this->m(), n = this->n(), n2b = this->n2b();
    auto border = this->border();
    for ( idxint i=0; i<m+2*border; i++ ) {
      for ( idxint j=0; j<n+2*border; j++ ) {
        char c = ( j<n+2*border-1 ? ' ' : '\n' );
        std::cout << std::format("{:5.2}{}",out[ this->oindex(i,j) ],c);
      }
    }
    return;
  };

  template< typename real >
  void bordered_array_sycl<real>::set_bc( bool down,bool right, bool trace ) {

  };

 template< typename real >
  void bordered_array_sycl<real>::scale_interior
      ( const sparsealg::bordered_array_base<real>& _other, real factor ) {
  };

};

namespace sparsealg {
  template class bordered_array_sycl<float>;
  template class bordered_array_sycl<double>;
};

