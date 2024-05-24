#include <cstddef>
#include <string>
#include <tuple>

std::tuple<bool,std::size_t,std::size_t,int,int,
	   bool, // gpu
	   bool, // trace
	   bool  // view
	   > parse_options
    (int,char**,std::string);
