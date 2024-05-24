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
    queue cpu_selector(cpu_selector_v);
    queue gpu_selector(gpu_selector_v);
    queue q;
    if (gpu)
        q = gpu_selector;
    else
        q = cpu_selector;

    std::cout << "Device : " << q.get_device().get_info<info::device::name>() << "\n";
    std::cout << "Max Compute Units : " << q.get_device().get_info<info::device::max_compute_units>() << std::endl;

    std::vector<real> Mat_A(msize*nsize,10.0);
    std::vector<real> Mat_Stencil(msize*nsize,0.0);

    buffer<real,2> Buf_a(Mat_A.data(),range<2>(msize,nsize));
    buffer<real,2> Buf_b(Mat_Stencil.data(),range<2>(msize,nsize));
    buffer<real,1> Buf_Fn(&FNorm, range<1>(1));

    using myclock = std::chrono::steady_clock;
    auto start_time = myclock::now();

    for (int repeat = 0; repeat<itcount; ++repeat) {

      // Kernel to compute the 5pt stencil and simultaneously the L2Norm
      q.submit([&] (handler &h)
      {
	accessor D_a(Buf_a,h);
	accessor D_b(Buf_b,h);
	auto D_Fn = reduction(Buf_Fn, h, std::plus<real>());

	h.parallel_for(range<2>(N-2,M-2), D_Fn, [=](item<2> index, auto &sum){
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
	
      q.submit([&] (handler &h)
      {
	accessor D_a(Buf_a,h);
	accessor D_b(Buf_b,h);
	accessor D_Fn(Buf_Fn,h);

	h.parallel_for(range<2>(N-2,M-2), [=](auto index){
	  auto row = index.get_id(0) + 1;
	  auto col = index.get_id(1) + 1;

	  D_a[row][col] = (D_b[row-1][col-1]/D_Fn[0]);
	});
      }).wait();

      end = rdtsc();
    
      elapsed_count[count] = (double)(end - start)/ClkPerSec;
      printf("TTC : %.12f\n",elapsed_count[count]);
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

    // for(int count=1;count<(atoi(argv[1]));count++)
    //   Average += elapsed_count[count];
   
    // std::cout << "\nTime to compute 5pt-Stencil + Power Method (Total) = " << Average << "\n"; 
    // Average = Average/(atoi(argv[1]) -1);
    // std::cout << "\nTime to compute (Avg over " << atoi(argv[1]) << " loops) = " << Average << "\n";

    return 0;
}
