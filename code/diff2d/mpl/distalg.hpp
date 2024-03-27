/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2023 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** distalg.hpp : headers for distributed bordered vector
 ****
 ****************************************************************/

#ifndef DISTALG_HPP
#define DISTALG_HPP

#define USE_SPAN

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <mpl/mpl.hpp>
#define BUFFER        std::pair<real*,std::shared_ptr<mpl::layout<real>>>
#define MKBUFFER std::make_pair<real*,std::shared_ptr<mpl::layout<real>>>

#include "../linalg.hpp"

namespace linalg {

  //codesnippet d2distarray
  template< typename real >
  class distributed_array {
  private:
    const mpl::cartesian_communicator& comm; 
    mpl::cartesian_communicator::dimensions dimensions; 
    //! orthogonal sizes of processor subdomains
    std::vector<size_t> proc_start_m,proc_start_n;
    //! coordinate of this processor
    mpl::cartesian_communicator::vector coord;
    std::map<char,int> neighbors;
    size_t m_global,n_global;
    bordered_array<real> data;
  //codesnippet end
  public:
    distributed_array
        ( const mpl::cartesian_communicator&,size_t m,size_t n,bool trace=false );
    std::pair<int,int> proc_grid_size();
    // std::pair<size_t,size_t> set_local_dimensions
    //     ( const mpl::cartesian_communicator&,int,size_t,size_t,bool=false );
    void set_neighbors( bool=false );

    BUFFER get_halo( char direction );
    BUFFER get_edge( char direction );
    void halo_exchange(bool=false);
    void central_difference_from( distributed_array<real>& other,bool=false );
    void copy_interior_from( distributed_array<real> other );
    void scale_interior( const distributed_array<real>&, real );
    real l2norm();
    void set( real value,bool trace=false) {
      data.set( value,trace );
    };
    void set_bc(bool trace=false) {
      data.set_bc( coord[0]==dimensions.size(0),coord[1]==dimensions.size(1), trace );
    };

    // logging
    void log_flops( float n );
    void log_bytes( float n );
    std::pair<float,float> log_report();

    void view( std::string="" );
  };
};

#endif
