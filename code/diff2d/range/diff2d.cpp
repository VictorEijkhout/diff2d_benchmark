/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2d.cpp : 2D diffusion in parallel through OpenMP
 ****
 ****************************************************************/

#include <chrono>
#include <format>
#include <iostream>
#include <memory>
using std::unique_ptr,std::make_unique;

#include "omp.h"
#include "../lib/range.hpp"
#include "../lib/options.hpp"

#include "tbb/tbb.h"

using real = float;

int main(int argc,char **argv) {

  auto [exit,msize,nsize,border,itcount,gpu,trace,view] = 
    parse_options(argc,argv,"Range algorithm version using execution policies");
  if (exit) return 0;

  const std::string prefix{"rng"};

  const int nthreads = 
    [] () {
    int nt;
#pragma omp parallel
#pragma omp single
    nt = omp_get_num_threads();
    return nt; }();
  tbb::global_control(tbb::global_control::max_allowed_parallelism, nthreads);

  std::cout << std::format("Runtype: {}\n",prefix);
  std::cout << std::format("Threads: {:>3}\n",nthreads);
  std::cout << std::format("Vector size: {} x {}\n",msize,nsize);

  auto X = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_range<real>>(msize,nsize,border) );
  auto Y = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_range<real>>(msize,nsize,border) );

  const int procno{0};
#include "../main.cpp"

  return 0;
}
