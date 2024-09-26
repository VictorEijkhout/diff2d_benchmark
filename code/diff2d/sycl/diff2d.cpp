/****************************************************************
 *  ****
 *  **** This file belongs with the course
 *  **** Parallel Programming in MPI and OpenMP
 *  **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 *  ****
 *  **** diff2d.cpp : 2D diffusion in parallel through SYCL
 *  **** based on code by Yojan Chitkara
 *  ****
 *  ****************************************************************/

#include <chrono>
#include <iostream>
using std::cout;
#include <format>
#include <memory>
using std::unique_ptr,std::make_unique;

#include "../lib/sycl.hpp"
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
  // codesnippet syclqueue
  sycl::queue q =
    [=] () -> sycl::queue {
      if (gpu) {
        cout << "Selecting device GPU\n";
        return queue(gpu_selector_v);
      } else {
        cout << "Selecting host CPU\n";
        return queue(cpu_selector_v);
      }
    }();
   //codesnippet syclqueue

  std::cout << std::format("Runtype: {}\n",prefix);
  std::cout << std::format("Threads: 1 obviously\n");
  std::cout << std::format("Vector size: {} x {}\n",msize,nsize);

  //codesnippet d2duniquesycl
  auto X = unique_ptr<linalg::bordered_array_base<real>>
      ( make_unique<linalg::bordered_array_sycl<real>>(msize,nsize,border,q) );
  //codesnippet end
  auto Y = unique_ptr<linalg::bordered_array_base<real>>
      ( make_unique<linalg::bordered_array_sycl<real>>(msize,nsize,border,q) );
  
  #include "../main.cpp"

  return 0;
}  
