/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
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
#define BUFFER        std::pair<real*,std::unique_ptr<mpl::layout<real>>>
#define MKBUFFER std::make_pair

#include "base.hpp"
#include "seq.hpp"

namespace linalg {

  //codesnippet d2distarray
  template< typename real >
  class distributed_array {
  private:
    size_t m_global,n_global;
    std::unique_ptr<bordered_array_base<real>> subdomain{nullptr};
    //codesnippet end

    //codesnippet d2distarraycomm    
    const mpl::cartesian_communicator& comm; 
    const int rank;
    mpl::cartesian_communicator::dimensions dimensions; 
    //! coordinate of this processor
    mpl::cartesian_communicator::vector coord;
    //codesnippet end    

    //codesnippet d2distarraydata
    // orthogonal sizes of processor subdomains
    std::vector<idxint> proc_start_m,proc_start_n;
    std::map<char,int> neighbors;
    // temp array, just for the central difference routine. somewhat wasteful
    std::unique_ptr<bordered_array_base<real>> tmp{nullptr};
  //codesnippet end
  public:
    // constructor
    distributed_array
        ( const mpl::cartesian_communicator&,size_t m,size_t n,int border,
          bool trace=false );
    std::vector<idxint> segmentize(idxint m,int pm,bool trace=false);
    void set_neighbors( bool=false );
    size_t global_n2b() const { return n_global+2*subdomain->border(); };
    std::tuple<size_t,size_t,int> outer_sizes() const {
      return std::make_tuple(m_global,n_global,subdomain->border()); };

    // required functionality
    void halo_exchange( bordered_array_base<real>&, bool=false);
    void central_difference_from
        ( const distributed_array<real>& other,bool=false );
    void scale_interior
        ( const linalg::distributed_array<real>& other, real );
    real l2norm() ;
    void set_value( real value,bool trace=false)  {
      subdomain->set_value( value,trace ); };
    void set_bc(bool down,bool right, bool trace=false)  {
      subdomain->set_bc( coord[0]==dimensions.size(0),coord[1]==dimensions.size(1), trace ); };
    void view( std::string="" ) ;

    /*
     * MPI communication routines
     */
    BUFFER get_halo( const bordered_array_base<real>&,char direction );
    BUFFER get_edge( const bordered_array_base<real>&,char direction );

    // logging
    void log_flops( float n );
    void log_bytes( float n );
    std::pair<float,float> log_report();

  };
};

#endif
