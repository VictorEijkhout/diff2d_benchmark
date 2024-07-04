/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2e.cpp : 2D diffusion in parallel using my own cartesian iterator
 ****
 ****************************************************************/

#include <chrono>
#include <iostream>
#include <format>
#include <memory>
using std::unique_ptr,std::make_unique;

#include "omp.h"
#include "../lib/options.hpp"

#include "diy2e.hpp"
using real = float;

int main(int argc,char **argv) {

  auto [exit,msize,nsize,border,itcount,gpu,trace,view] =
    parse_options(argc,argv,"OpenMP over 1D index range with mdspan");
  if (exit) return 0;
  const std::string prefix{"diy"};
  const int procno{0};

#pragma omp parallel
#pragma omp single
  std::cout << std::format("Total cores available: {}\n",omp_get_num_procs());
  
  const int nthreads = 
    [] () {
    int nt;
#pragma omp parallel
#pragma omp single
    nt = omp_get_num_threads();
    return nt; }();

  std::cout << std::format("Runtype: {}\n",prefix);
  std::cout << std::format("Threads: {:>3}\n",nthreads);
  std::cout << std::format("Vector size: {} x {}\n",msize,nsize);

  auto X = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_diy2e<real>>(msize,nsize,border) );
  auto Y = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_diy2e<real>>(msize,nsize,border) );

#include "../main.cpp"

  return 0;
}
