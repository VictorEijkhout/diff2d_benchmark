/****************************************************************
 *  ****
 *  **** This file belongs with the course
 *  **** Parallel Programming in MPI and OpenMP
 *  **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
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

//codesnippet d2dsyclindex
#define IINDEX( i,j,m,n,b ) ((i)+b)*(n+2*b) + (j)+b
//codesnippet end

namespace linalg {
  template< typename real >
  bordered_array_sycl<real>::bordered_array_sycl( idxint m,idxint n,int border, queue q )
    : bordered_array_base<real>(m,n,border) {
    auto out = this->data();
    //[m,n,border] = this->outer_sizes();

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
      ( const linalg::bordered_array_base<real>& _other,bool trace ) {

    const auto& other = dynamic_cast<const linalg::bordered_array_sycl<real>&>(_other);

    auto out = this->data();
    auto in = other.data();
    auto [m,n,b,m2b,n2b] = this->inner_sizes();
    auto q = this->q;
    buffer<real,2> Buf_a(out, sycl::range<2>(m2b,n2b));
    buffer<real,2> Buf_b(in, sycl::range<2>(m,n));
    //codesnippet d2d5ptsycl
    q.submit([&] (handler &h)
    {
      accessor D_a(Buf_a,h,read_only);
      accessor D_b(Buf_b,h,write_only);

      h.parallel_for(range<2>(other.m() - 2,other.n() - 2), [=](item<2> index){
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;

        real stencil_value =
          4*D_a[row][col]
          - D_a[row-1][col] - D_a[row+1][col]
          - D_a[row][col-1] - D_a[row][col+1];
        D_b[row-1][col-1] = stencil_value;
      });
    }).wait();
    //codesnippet end

  };

  template < typename real >
  void bordered_array_sycl<real>::set( real value, bool trace ) {
    auto out = this->data();
    //auto [m,n,b,m2b,n2b] = this->inner_sizes();
    auto m = std::get<0>(this->inner_sizes());
    auto n = std::get<1>(this->inner_sizes());

    auto q = this->q;
    //codesnippet syclbufaccess
    buffer<real,2> Buf_a(out, sycl::range<2>(m,n));
    q.submit([&](sycl::handler &h) {
      sycl::accessor D_a(Buf_a, h, sycl::write_only);
    //codesnippet end

      h.parallel_for(sycl::range<2>(m-2, n-2), [=](auto index) {
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;
        D_a[row][col] = value;
      });
    }).wait();
  };

  //codesnippet d2dnormsycl
  template < typename real >
  real bordered_array_sycl<real>::l2norm() {
    real norm = 0.0;
    auto out = this->data();
    auto m = std::get<0>(this->inner_sizes());
    auto n = std::get<0>(this->inner_sizes());
    buffer<real,2> Buf_a(out, sycl::range<2>(m,n));
    buffer<real,1> Buf_n(&norm, range<1>(1));

    auto q = this->q;
    q.submit([&](handler &h) {
      accessor D_a(Buf_a, h, read_only);
      auto D_Fn = reduction(Buf_n, h, std::plus<real>());

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
  };

  template< typename real >
  void bordered_array_sycl<real>::set_bc( bool down,bool right, bool trace ) {

  };

 template< typename real >
  void bordered_array_sycl<real>::scale_interior
      ( const linalg::bordered_array_base<real>& _other, real factor ) {
  };

};

namespace linalg {
  template class bordered_array_sycl<float>;
  template class bordered_array_sycl<double>;
};

