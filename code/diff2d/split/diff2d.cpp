/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2d.cpp : 2D diffusion in parallel through OpenMP
 ****
 ****************************************************************/

#define USE_OMP

#include <chrono>

#include <algorithm>
using std::for_each;

#include <cassert>

#include <iostream>
using std::cout;
#include <format>
using std::format;

#include "cxxopts.hpp"

#include "omp.h"
#include "../linalg.hpp"

using real = float;

int main(int argc,char **argv) {

#include "../options.cpp"
  if (result.count("Help")) {
    std::cout << "================ This type of run:\n";
    std::cout << "Very inefficient 1D range split into 2D i,j\n";
    std::cout << "indexing into 2D mdspan\n";
    std::cout << "parallel by OpenMP\n";
    return 0;
  }

  const int nthreads = 
    [] () {
    int nt;
#pragma omp parallel
#pragma omp single
    nt = omp_get_num_threads();
    return nt; }();
  const std::string prefix{"span"};

  cout << format("Threads: {:>3}\n",nthreads);
  if (trace)
    cout << format("Vector size: {} x {}\n",msize,nsize);
  linalg::bordered_array<real> X(msize,nsize,border),Y(msize,nsize,border);
  X.set_collapse(collapse); Y.set_collapse(collapse);

  // force tracing calls to happen
  const int procno{0};
#include "../main.cpp"

  return 0;
}
