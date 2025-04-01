/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2d.cpp : 2D diffusion in parallel through MPL
 ****
 ****************************************************************/

#include <chrono>
#include <iostream>
#include <format>
#include <memory>
using std::unique_ptr,std::make_unique;

#include <cassert>

#include <mpl/mpl.hpp>
#include "../lib/dist.hpp"
#include "../lib/options.hpp"

using real = float;

int main(int argc,char **argv) {

  const mpl::communicator &comm_world=mpl::environment::comm_world();
  auto procno = comm_world.rank();
  auto nprocs = comm_world.size();

  auto [exit,msize,nsize,border,itcount,gpu,trace,view] =
    parse_options(argc,argv,"OpenMP version using 2D loop");
  comm_world.bcast(0,exit);
  if (exit) return 0;
  const std::string prefix{"mpl"};
  comm_world.bcast(0,msize);
  comm_world.bcast(0,nsize);
  comm_world.bcast(0,border);
  comm_world.bcast(0,itcount);
  comm_world.bcast(0,trace);
  comm_world.bcast(0,view);

  if (procno==0) {
    std::cout << std::format("Runtype: {}\n",prefix);
    std::cout << std::format("Processes: {:>3}\n",nprocs);
    std::cout << std::format("Vector size: {} x {}\n",msize,nsize);
  }
  
  //codesnippet d2ddimscreate
  mpl::cartesian_communicator::dimensions brick(2);
  brick = mpl::dims_create(nprocs,brick);
  mpl::cartesian_communicator cart_comm( comm_world,brick );
  //codesnippet end

  if ( trace and procno==0 )
    std::cout << std::format("Process grid: {}x{}\n",brick.size(0),brick.size(1));
  if ( trace and procno==0 )
    std::cout << std::format("Setting up domain of {} x {}\n",msize,nsize);

  auto X = // unique_ptr<sparsealg::bordered_array_base<real>>
    make_unique<sparsealg::distributed_array<real>>(cart_comm,msize,nsize,border,trace);
  auto Y = // unique_ptr<sparsealg::bordered_array_base<real>>
    make_unique<sparsealg::distributed_array<real>>(cart_comm,msize,nsize,border,trace);

  try {
#include "../main.cpp"
  } catch ( const char* c ) {
    std::cout << std::format("Exception: {}\n",c);
  } catch ( const int i ) {
    std::cout << std::format("Exception: {}\n",i);
  } catch ( ... ) {
    std::cout << std::format("Exception of unknown type\n");
  }

  return 0;
}
