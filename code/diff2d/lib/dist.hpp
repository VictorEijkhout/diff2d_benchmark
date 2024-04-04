/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** dist.hpp : headers for distributed bordered vector
 ****
 ****************************************************************/

#ifndef DISTALG_HPP
#define DISTALG_HPP

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <mpl/mpl.hpp>
#define BUFFER        std::pair<real*,std::shared_ptr<mpl::layout<real>>>
#define MKBUFFER std::make_pair<real*,std::shared_ptr<mpl::layout<real>>>

#include "base.hpp"
#include "seq.hpp"

namespace linalg {

  //codesnippet d2distarray
  template< typename real >
  class distributed_array { // : public bordered_array_base<real> {
  private:
    const mpl::cartesian_communicator& comm; 
    mpl::cartesian_communicator::dimensions dimensions; 
    //! orthogonal sizes of processor subdomains
    std::vector<size_t> proc_start_m,proc_start_n;
    //! coordinate of this processor
    mpl::cartesian_communicator::vector coord;
    std::map<char,int> neighbors;
    size_t m_global,n_global;
    std::unique_ptr<bordered_array_base<real>> data{nullptr},tmp{nullptr};
  //codesnippet end
  public:
    // constructor
    distributed_array
        ( const mpl::cartesian_communicator&,size_t m,size_t n,bool trace=false );
    void set_neighbors( bool=false );
    size_t global_n2b() const { return n_global+2*data->border(); };

    // required functionality
    void central_difference_from
        ( const linalg::bordered_array_base<real>&,bool=false );
    void scale_interior( const linalg::bordered_array_base<real>&, real );
    real l2norm() ;
    void set( real value,bool trace=false)  {
      data->set( value,trace ); };
    void set_bc(bool down,bool right, bool trace=false)  {
      data->set_bc( coord[0]==dimensions.size(0),coord[1]==dimensions.size(1), trace ); };
    void view( std::string="" ) ;

    /*
     * MPI communication routines
     */
    BUFFER get_halo( char direction );
    BUFFER get_edge( char direction );
    void halo_exchange( bordered_array_base<real>&,bool=false);
    void central_difference_from( const distributed_array<real>& other,bool=false );
    void copy_interior_from( const distributed_array<real> other );
    void scale_interior( const distributed_array<real>&, real );

    // logging
    void log_flops( float n );
    void log_bytes( float n );
    std::pair<float,float> log_report();

  };
};

#endif
