/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2d.cpp : 2D diffusion in parallel through SYCL
 **** based on code by Yojan Chitkara
 ****
 ****************************************************************/

#include <chrono>
#include <iostream>
using std::cout;
#include <format>

#include <sycl/sycl.hpp>
using namespace sycl;

#include "../lib/options.hpp"
using real = float;

int main(int argc,char *argv[])
{

  auto [exit,msize,nsize,border,itcount,gpu,trace,view] =
    parse_options(argc,argv,"SYCL version using 2D loop");
  if (exit) return 0;
  const std::string prefix{"sycl"};
  const int procno{0};

  /*
   * SYCL setup
   */
    queue q =
      [=] () -> queue {
	if (gpu) {
	  cout << "Selecting device GPU\n";
	  return queue(gpu_selector_v);
	} else {
	  cout << "Selecting host CPU\n";
	  return queue(cpu_selector_v);
	}
      }();

    std::cout << "Device : " << q.get_device().get_info<info::device::name>() << "\n";
    std::cout << "Max Compute Units : " << q.get_device().get_info<info::device::max_compute_units>() << std::endl;

    std::vector<real> Mat_A(msize*nsize,10.0);
    std::vector<real> Mat_Stencil(msize*nsize,0.0);
    real FNorm;

    buffer<real,2> Buf_a(Mat_A.data(),range<2>(msize,nsize));
    buffer<real,2> Buf_b(Mat_Stencil.data(),range<2>(msize,nsize));
    buffer<real,1> Buf_Fn(&FNorm, range<1>(1));

    q.submit([&] (handler &h) {
      accessor D_a(Buf_a,h);
      accessor D_b(Buf_b,h);
      accessor D_Fn(Buf_Fn,h);

      h.parallel_for(range<2>(msize-2,nsize-2), [=](auto index){
	auto row = index.get_id(0) + 1;
	auto col = index.get_id(1) + 1;

	D_a[row][col] = 1.;
      });
    }).wait();

    using myclock = std::chrono::steady_clock;
    auto start_time = myclock::now();

    for (int it = 0; it<itcount; ++it) {

      // Kernel to compute the 5pt stencil and simultaneously the L2Norm
      q.submit([&] (handler &h)
      {
	accessor D_a(Buf_a,h);
	accessor D_b(Buf_b,h);
	auto D_Fn = reduction(Buf_Fn, h, std::plus<real>());

	h.parallel_for(range<2>(msize-2,nsize-2), D_Fn, [=](item<2> index, auto &sum){
	  auto row = index.get_id(0) + 1;
	  auto col = index.get_id(1) + 1;

	  real stencil_value =
	    4*D_a[row][col]
	    - D_a[row-1][col] - D_a[row+1][col]
	    - D_a[row][col-1] - D_a[row][col+1];
	  D_b[row-1][col-1] = stencil_value;
	  sum += (stencil_value * stencil_value);
	});
      }).wait();

      FNorm = std::sqrt(FNorm);
      if (trace)
	std::cout << std::format("[{:>2}] y norm: {}\n",it,FNorm);
	
      q.submit([&] (handler &h)
      {
	accessor D_a(Buf_a,h);
	accessor D_b(Buf_b,h);
	accessor D_Fn(Buf_Fn,h);

	h.parallel_for(range<2>(msize-2,nsize-2), [=](auto index){
	  auto row = index.get_id(0) + 1;
	  auto col = index.get_id(1) + 1;

	  D_a[row][col] = (D_b[row-1][col-1]/D_Fn[0]);
	});
      }).wait();

    }

    /*
     * Time reporting
     */
    auto duration = myclock::now()-start_time;
    auto millisec_duration = 
      std::chrono::duration_cast<std::chrono::microseconds>(duration)/1000;
    auto msec = millisec_duration.count();
    if ( procno==0 )
      std::cout << std::format("Time: {:>6} msec\n",msec);

    return 0;
}
