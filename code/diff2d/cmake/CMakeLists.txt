cmake_minimum_required( VERSION 3.20 )
if( NOT DEFINED VARIANT )
  message( FATAL_ERROR "\nYou need to set the VARIANT variable\n" )
endif()

project( ${VARIANT} VERSION 1.0 )
message( STATUS "\n====   Building variant: ${PROJECT_NAME}" )

set( CMAKE_BUILD_TYPE RelWithDebInfo )
add_executable( ${PROJECT_NAME} )
target_compile_features( ${PROJECT_NAME} PRIVATE cxx_std_23 )
## set_target_properties(${PROJECT_NAME}_c PROPERTIES CXX_STANDARD 23)
target_sources( ${PROJECT_NAME} PRIVATE diff2d.cpp )
target_sources( ${PROJECT_NAME} INTERFACE ../main.cpp )

##
## parallelism
##
message( STATUS "\n====   setting backend" )
# what variant are we doing?
set( omp_backends "seq" "oned" "clps" "iota" "diy2e" "kokkos" "span" )
list( FIND omp_backends "${VARIANT}" index )
# add OpenMP 
if( ${index} GREATER 0 )
    message( STATUS "\n====   adding OpenMP backend for index ${index}" )
    find_package(OpenMP)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX  )
elseif( ${VARIANT} STREQUAL "sycl" )
  message( STATUS "\n====   using Sycl backend" )
  ## find_package( IntelDPCPP REQUIRED )
  find_package( IntelSYCL REQUIRED )
  add_sycl_to_target(TARGET ${PROJECT_NAME} SOURCES diff2d.cpp )
endif()


include( FetchContent )

##
## kokkos handling
##

if( ${PROJECT_NAME} STREQUAL "kokkos" )
  find_package( Kokkos QUIET )
  if( Kokkos_FOUND )
    message( STATUS "Thanks for installing kokkos yourself" )
  else()
    message( STATUS "Downloading kokkos for you" )
    FetchContent_Declare(
      Kokkos
      GIT_REPOSITORY https://github.com/kokkos/kokkos.git
      )
    FetchContent_MakeAvailable( Kokkos )
  endif()
  target_link_libraries( ${PROJECT_NAME} PUBLIC Kokkos::kokkos )
endif()

##
## options handling
##

find_package( cxxopts QUIET )
if( cxxopts_FOUND )
  message( STATUS "Thanks for installing cxxopts yourself" )
else()
  message( STATUS "Downloading cxxopts for you" )
  FetchContent_Declare(
    cxxopts GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git )
  FetchContent_MakeAvailable( cxxopts )
endif()
target_include_directories( ${PROJECT_NAME} PUBLIC cxxopts::cxxopts )
target_link_libraries( ${PROJECT_NAME} PUBLIC cxxopts::cxxopts )
target_sources( ${PROJECT_NAME} PRIVATE ../lib/options.cpp )

##
## span handling
##
find_package( mdspan QUIET )
if( mdspan_FOUND )
  message( STATUS "Thanks for installing mdspan yourself" )
else()
  message( STATUS "Downloading mdspan for you" )
  FetchContent_Declare(
    mdspan
    GIT_REPOSITORY https://github.com/kokkos/mdspan.git
    GIT_TAG stable
    )
  FetchContent_MakeAvailable( mdspan )
endif()
target_link_libraries( ${PROJECT_NAME} PUBLIC std::mdspan )

target_sources( ${PROJECT_NAME} PRIVATE ../lib/${VARIANT}.cpp ../lib/base.cpp )
target_include_directories( ${PROJECT_NAME} PRIVATE ../lib )

install( TARGETS ${PROJECT_NAME} DESTINATION . )
