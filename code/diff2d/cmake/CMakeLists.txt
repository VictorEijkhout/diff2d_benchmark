cmake_minimum_required( VERSION 3.20 )
if( NOT DEFINED VARIANT )
  message( FATAL_ERROR "\nYou need to set the VARIANT variable\n" )
endif()

project( ${VARIANT} VERSION 1.0 )
message( STATUS "\n====   Building variant: ${PROJECT_NAME}" )

add_executable( ${PROJECT_NAME} )
target_compile_features( ${PROJECT_NAME} PRIVATE cxx_std_23 )
target_sources( ${PROJECT_NAME} PRIVATE diff2d.cpp )
target_sources( ${PROJECT_NAME} INTERFACE ../main.cpp )

##
## parallelism
##
message( STATUS "\n====   setting backend" )
# what variant are we doing?
set( omp_backends "seq" "oned" "clps" "kokkos" "span" )
list( FIND omp_backends "${VARIANT}" index )
# add OpenMP 
if( ${index} GREATER 0 )
    message( STATUS "\n====   adding OpenMP backend for index ${index}" )
    find_package(OpenMP)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX  )
endif()
# if( ${index} EQUAL 3 )
#   message( STATUS "   also adding Kokkos" )
#   # Create imported target std::mdspan
#   target_link_libraries( ${PROJECT_NAME} std::mdspan)
#   add_library(std::mdspan INTERFACE IMPORTED)
#   set_target_properties(std::mdspan PROPERTIES
#     INTERFACE_INCLUDE_DIRECTORIES ${MDSPAN_INC} )
# endif()

##
## options handling
##

message( STATUS "\n====   options handling" )
# include( FetchContent )
# FetchContent_Declare(
#     cxxopts GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git )
# FetchContent_MakeAvailable( cxxopts )
# target_include_directories( ${PROJECT_NAME} PUBLIC cxxopts::cxxopts )
# target_link_libraries( ${PROJECT_NAME} PUBLIC cxxopts::cxxopts )

find_package( cxxopts REQUIRED )
target_link_libraries( ${PROJECT_NAME} PUBLIC cxxopts::cxxopts )
target_sources( ${PROJECT_NAME} PRIVATE ../lib/options.cpp )

##
## span handling
##
message( STATUS "\n====   span handling" )
find_package( mdspan REQUIRED )
target_link_libraries( ${PROJECT_NAME} PUBLIC std::mdspan )
target_sources( ${PROJECT_NAME} PRIVATE ../lib/${VARIANT}.cpp ../lib/base.cpp )
target_include_directories( ${PROJECT_NAME} PRIVATE ../lib )

install( TARGETS ${PROJECT_NAME} DESTINATION . )
