/****************************************************************
 *  ****
 *  **** This file belongs with the course
 *  **** Parallel Programming in MPI and OpenMP
 *  **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 *  ****
 *  **** sycl.cpp : 2D diffusion in parallel through SYCL 
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
    const auto [m,n,b] = this->outer_sizes();

    auto q = this->q;
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
    //codesnippet d2d5ptsycl
    q.submit([&] (handler &h)
    {
      accessor D_a(Buf_a,h,read_only);
      accessor D_b(Buf_b,h,write_only);

      h.parallel_for(range<2>(other.m - 2,other.n - 2) [=](item<2> index){
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;

        real stencil_value =
          4*D_a[row][col]
          - D_a[row-1][col] - D_a[row+1][col]
          - D_a[row][col-1] - D_a[row][col+1];
        D_b[row-1][col-1] = stencil_value;
      });
    }).wait();    
  };

  template < typename real >
  void bordered_array_sycl<real>::set( real value, bool trace ) {
    auto &A = this;
    buffer<real,2> Buf_a(A.internal_data().data(), sycl::range<2>(A.m(),A.n()));
    auto q = this->q;

    q.submit([&](sycl::handler &h) {
      sycl::accessor D_a(Buf_a, h, sycl::write_only);

      h.parallel_for(sycl::range<2>(A.m()-2, A.n()-2), [=](auto index) {
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;
        D_a[row][col] = value;
      });
    }).wait();
  };

  template < typename real >
  real bordered_array_sycl<real>::l2norm() {
    real norm = 0.0;
    auto &A = this;
    buffer<real,2> Buf_a(A.internal_data().data(), sycl::range<2>(A.m(),A.n()));
    buffer<real,1> Buf_n(&norm, range<1>(1));

    auto q = this->q;
    q.submit([&](handler &h) {
      accessor D_a(Buf_a, h, read_only);
      auto D_f = sycl::reduction(&Buf_n, std::plus<real>());

      h.parallel_for(range<2>(A.m() - 2, A.n() - 2), D_f, [=](item<2> index, auto &sum) {
        auto row = index.get_id(0) + 1;
        auto col = index.get_id(1) + 1;

        sum += (D_a[row-1][col-1] * D_a[row-1][col-1]);
      });
    }).wait();    

    norm = std::sqrt(norm);
    return norm;
  };

  template < typename real >
  void bordered_array_sycl<real>::view( string caption ) {
    if (caption!="")
      cout << format("{}:\n",caption);

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
};

namespace linalg {
  template class bordered_array_sycl<float>;
  template vlass bordered_array_sycl<double>
};
