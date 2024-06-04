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
set( omp_backends "seq" "oned" "clps" "kokkos" "span" )
list( FIND omp_backends "${VARIANT}" index )
if( ${index} GREATER 1 )
else()
    message( STATUS "\n====   adding OpenMP backend" )
    find_package(OpenMP)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX )
endif()

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

find_package( PkgConfig REQUIRED )
pkg_check_modules( OPTS REQUIRED cxxopts )
target_include_directories( ${PROJECT_NAME} PRIVATE ${OPTS_INCLUDE_DIRS} )
target_sources( ${PROJECT_NAME} PRIVATE ../lib/options.cpp )

##
## span handling
##
message( STATUS "\n====   span handling" )
target_sources( ${PROJECT_NAME} PRIVATE ../lib/${VARIANT}.cpp ../lib/base.cpp )
target_include_directories( ${PROJECT_NAME} PRIVATE ../lib )
find_package( mdspan REQUIRED )
## target_link_libraries( ${PROJECT_NAME} PRIVATE mdspan::mdspan )
## target_include_directories( ${PROJECT_NAME} PRIVATE mdspan::mdspan )
target_include_directories( ${PROJECT_NAME} PRIVATE ${MDSPAN_INC} )

install( TARGETS ${PROJECT_NAME} DESTINATION . )
