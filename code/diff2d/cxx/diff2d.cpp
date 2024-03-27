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
#include <iostream>
#include <format>
#include <memory>
using std::unique_ptr,std::make_unique;

#include "../lib/seq.hpp"
#include "../lib/options.hpp"

using real = float;

int main(int argc,char **argv) {

  auto [exit,msize,nsize,border,itcount,trace,view] =
    parse_options(argc,argv,"OpenMP version using 2D loop");
  if (exit) return 0;
  const std::string prefix{"seq"};
  const int procno{0};

  std::cout << std::format("Runtype: {}\n",prefix);
  std::cout << std::format("Threads: 1 obviously\n");
  std::cout << std::format("Vector size: {} x {}\n",msize,nsize);

  auto X = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_seq<real>>(msize,nsize,border) );
  auto Y = unique_ptr<linalg::bordered_array_base<real>>
    ( make_unique<linalg::bordered_array_seq<real>>(msize,nsize,border) );

#include "../main.cpp"

  return 0;
}
