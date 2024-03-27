#include <cstddef>
#include <string>
#include <tuple>

std::tuple<bool,std::size_t,std::size_t,int,int,bool,bool> parse_options
    (int,char**,std::string);
