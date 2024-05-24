/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** options.cpp : options handling for diff2d
 ****
 ****************************************************************/

#include <cstddef>
using std::size_t;
#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "options.hpp"

std::tuple<bool,size_t,size_t,int,int,bool,bool> parse_options
    (int argc,char **argv,std::string helpmsg) {
  cxxopts::Options options
    ("diff2d",
     "================================\nPrime numbers\n================================");
  options.add_options()
    ("h,help","usage information")
    ("H,Help","elaborate information about this type of run")
    ("m,msize","size in 1st direction",cxxopts::value<size_t>()->default_value("10"))
    ("n,nsize","size in 2nd direction",cxxopts::value<size_t>()->default_value("10"))
    ("b,border","border around omega",cxxopts::value<int>()->default_value("1"))
    ("i,itcount","iteration count",cxxopts::value<int>()->default_value("10"))
    ("p,parallelism","number of threads/processes",cxxopts::value<int>()->default_value("4"))
    ("g,gpu","use GPU",cxxopts::value<bool>()->default_value("false"))
    ("t,trace","view arrays",cxxopts::value<bool>()->default_value("false"))
    ("v,view","view arrays",cxxopts::value<bool>()->default_value("false"))
    ;

  auto result = options.parse(argc,argv);
  if (result.count("help")) {
    std::cout << options.help() << '\n';
    return std::make_tuple(true,0,0,0,0,false,false);
  }

  if (result.count("Help")) {
    std::cout << "================ This type of run:\n";
    std::cout << helpmsg << '\n';
    return std::make_tuple(true,0,0,0,0,false,false);
  }

  size_t msize    = result["m"].as<size_t>();
  size_t nsize    = result["n"].as<size_t>();
  int border   = result["b"].as<int>();
  int itcount  = result["i"].as<int>();

  bool gpu     = result["g"].as<bool>();
  bool trace   = result["t"].as<bool>();
  bool view    = result["v"].as<bool>();
  bool help_msg = result.count("Help")>0;

  return std::make_tuple(false,msize,nsize,border,itcount,gpu,trace,view);
}
