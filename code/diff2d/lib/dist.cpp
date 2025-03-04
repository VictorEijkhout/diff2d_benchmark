/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** dist.cpp : distributed vector routines
 ****
 ****************************************************************/

#include <memory>
using std::unique_ptr,std::make_unique;
#include <string>
using std::string;
#include <tuple>
using std::make_pair;
#include <vector>
using std::vector;

#include <iostream>
using std::cout;
#include <format>
using std::format;

#include <ranges>
namespace rng = std::ranges;

#include "dist.hpp"

#include <memory>
using std::unique_ptr,std::make_unique;

namespace linalg {

  /*
   * Constructor
   */
  //codesnippet d2ddistributed1
  template< typename real >
  distributed_array<real>::distributed_array
      ( const mpl::cartesian_communicator& comm,size_t m,size_t n,int border,
        bool trace )
        : comm(comm)
        , procrank( comm.rank() )
        , m_global(m),n_global(n)
  {
    // compute local sizes...
    //codesnippet end
    /*
     * set up processor grid
     */
    //codesnippet d2ddistributed2
    // global processor grid size
    dimensions = comm.get_dimensions();
    auto
      pmsize = dimensions.size(0),
      pnsize = dimensions.size(1);
    //codesnippet end

    /*
     * Start and end for each processor
     */
    proc_start_m = segmentize(m,pmsize); 
    proc_start_n = segmentize(n,pnsize);
    /*
     * This process in particular
     */
    //codesnippet d2ddistributed3
    coord = comm.coordinates(procrank);
    auto pm = coord[0], pn = coord[1];
    auto
      sm = proc_start_m[pm+1]-proc_start_m[pm],
      sn = proc_start_n[pn+1]-proc_start_n[pn];
    //codesnippet end
    //codesnippet d2ddistributed4
    subdomain = unique_ptr<bordered_array_base<real>>
      ( make_unique<bordered_array_seq<real>>(sm,sn,border) );
    //codesnippet end
    tmp = unique_ptr<bordered_array_base<real>>
      ( make_unique<bordered_array_seq<real>>(sm,sn,border) );
    set_neighbors(trace);
  };

  
  template< typename real >
  std::vector<idxint> distributed_array<real>::segmentize
      (idxint size,int psize,bool trace) {
    vector<idxint> segments(psize+1);
    for ( int pi=0; pi<=psize; pi++ )
      segments.at(pi) = pi*size/psize;
    if (trace and procrank==0) 
      rng::for_each( segments,
                     [init=true] ( auto p ) mutable {
                       if (init) { init=false; cout << format("pstarts: {}",p); } 
                       else cout << format(", {}",p); } ); cout << format("\n");
    return segments;
  };

  /*! Auxiliary for constructor.
   * We use matrix-numbering: down and right
   */
  //codesnippet d2dnsew
  template< typename real >
  void distributed_array<real>::set_neighbors( bool trace ) {
    mpl::shift_ranks shifted = comm.shift
      ( /* dim: */ 0,/* increasing: */ 1 );
    neighbors['S'] = 
      [shifted] () { auto v = shifted.destination; 
                     return (v<0 ? mpl::proc_null : v); }();
    neighbors['N'] = 
      [shifted] () { auto v = shifted.source; 
                     return (v<0 ? mpl::proc_null : v); }();
    shifted = comm.shift
      ( /* dim: */ 1,/* increasing: */ 1 );
    neighbors['E'] = 
      [shifted] () { auto v = shifted.destination;
                     return (v<0 ? mpl::proc_null : v); }();
    neighbors['W'] = 
      [shifted] () { auto v = shifted.source; 
                     return (v<0 ? mpl::proc_null : v); }();
    //codesnippet end
    if (trace) { 
      cout << format("[{}] ",comm.rank());
      for_each ( neighbors.begin(),neighbors.end(),
                 [init=true] ( auto n ) mutable {
                   if (init) init=false; else cout << format(", ");
                   auto [d,r] = n; cout << format("{}={}",d,r);
                 } );
      cout << format("\n");
    }
  };

  /*! Derive pointer and layout for the inner edge send buffer, 
   * given the target direction */
  template< typename real >
  BUFFER distributed_array<real>::get_edge( const bordered_array_base<real> &data,char direction ) {
    const auto [m,n,b,m2b,n2b] = data.inner_sizes();
    const auto data2d = data.data2d();
    real* ptr; 
    if ( direction=='N' )
      ptr = &( data2d[b,b] );
    else if ( direction=='W' )
      ptr = &( data2d[b,b] );
    else if ( direction=='E' )
      ptr = &( data2d[b,n+b-1] );
    else if ( direction=='S' )
      ptr = &( data2d[m+b-1,b] );
    unique_ptr<mpl::layout<real>> layout;
    if (b>1)
      throw( "can not make layout for b>1" );
    if ( direction=='N' or direction=='S' )
      //codesnippet d2dlythor
      layout = unique_ptr<mpl::layout<real>>
        ( make_unique<mpl::contiguous_layout<real>>(n) );
      //codesnippet end
    else
      //codesnippet d2dlytver
      layout = unique_ptr<mpl::layout<real>>
        ( make_unique<mpl::strided_vector_layout<real>>(m,1,n+2*b) );
      //codesnippet end
    return make_pair(ptr,layout);      
  };

  /*! Derive pointer and layout for the halo receive buffer,
   *  given the source direction */
  template< typename real >
  BUFFER distributed_array<real>::get_halo( const bordered_array_base<real> &data,char direction ) {
    // implicit assumption that `data' is compatible with `this'
    const auto [m,n,b,m2b,n2b] = data.inner_sizes();
    const auto data2d = data.data2d();
    real* ptr; 
    if ( direction=='N' )
      ptr = &( data2d[0,b] );
    else if ( direction=='W' )
      ptr = &( data2d[b,0] );
    else if ( direction=='E' )
      ptr = &( data2d[b,n+b] );
    else if ( direction=='S' )
      ptr = &( data2d[m+b,b] );
    unique_ptr<mpl::layout<real>> layout;
    if (b>1)
      throw( "can not make layout for b>1" );
    if ( direction=='N' or direction=='S' )
      layout = unique_ptr<mpl::layout<real>>
        ( make_unique<mpl::contiguous_layout<real>>(n) );
    else
      layout = unique_ptr<mpl::layout<real>>
        ( make_unique<mpl::strided_vector_layout<real>>(m,1,n+2*b) );
    return MKBUFFER(ptr,layout);      
  };

  /*! Halo exchange by loopoing over 4 directions.
   * This is almost a class method: it operates completely on a bordered array
   * passed in as parameter.
   */
  //codesnippet d2dexchange
  template< typename real >
  void distributed_array<real>::halo_exchange
      ( bordered_array_base<real>& data, bool trace ) {
    for ( auto [f,t] : { make_pair('W','E'),{'E','W'},{'N','S'},{'S','N'} } ) {
      // inner edge to send in the `to' direction
      auto [snd_ptr,snd_layout] = get_edge(data,t);
      // halo region to receive from the `from' direction
      auto [rcv_ptr,rcv_layout] = get_halo(data,f);
      // neighbor ranks
      auto from = neighbors[f], to = neighbors[t];
      mpl::tag_t t0{0};
      comm.sendrecv
        ( snd_ptr,*snd_layout,to,  t0,
          rcv_ptr,*rcv_layout,from,t0);
    }
  };
  //codesnippet end

  //! Other is by non-const reference because we do halo exchange on it.
  template< typename real >
  void distributed_array<real>::central_difference_from
  ( const linalg::distributed_array<real>& other,bool trace ) {
    // make distributed vector that has internal data of `other'
    tmp->scale_interior( *(other.subdomain),static_cast<real>(1) );
    halo_exchange( *tmp,trace );
    if (trace)
      tmp->view("halo exchanged");
    subdomain->central_difference_from( *tmp,trace );
  };

  //! Norm computation
  //codesnippet d2dnormmpl
  template< typename real >
  real distributed_array<real>::l2norm() { 
    real local_norm = subdomain->l2norm();
    local_norm *= local_norm;
    real global_norm;
    comm.allreduce( mpl::plus<real>(),local_norm,global_norm);
    return std::sqrt(global_norm);
  };
  //codesnippet end

  //! Scale the local array
  //codesnippet d2dscalempl
  template< typename real >
  void distributed_array<real>::scale_interior
      ( const linalg::distributed_array<real>& other, real factor ) {
    subdomain->scale_interior( *(other.subdomain),factor);
  };
  //codesnippet end

  //! Gather an array for printing
  template< typename real >
  void distributed_array<real>::view( string caption ) {
    size_t global_size{0}; auto myrank = comm.rank();
    // size without borders
    size_t local_size = subdomain->inner_size();
    // send data is strided vector because of border
    real* first_data_point =  // should this be an mdspan?
      subdomain->data() + 2*subdomain->n2b() + 2*subdomain->border();
    auto proc_inner_layout =
      [ this ] ( auto procrank,auto n2b ) { 
        auto coord = this->comm.coordinates(procrank);
        auto pm = coord[0], pn = coord[1];
        auto msize = this->proc_start_m[pm+1]-this->proc_start_m[pm];
        auto bs    = this->proc_start_n[pn+1]-this->proc_start_n[pn];
        return mpl::strided_vector_layout<real>( msize, bs, n2b );
      };
    mpl::layout<real> send_layout = proc_inner_layout(comm.rank(),subdomain->n2b());
    subdomain->view( std::format("[{}] {}",myrank,caption) );
    return;
    if ( myrank==0 ) {
      comm.reduce( mpl::plus<size_t>(),0,local_size,global_size );
      vector<real> gathered_data(global_size);
      const auto nprocs = comm.size();
      mpl::layouts<real> recv_layouts(nprocs);
      mpl::displacements recv_displs(nprocs);
      for ( auto iproc=0; iproc<nprocs; iproc++ ) {
        auto coord = comm.coordinates(iproc);
        auto pm = coord[0], pn = coord[1];
        auto msize = proc_start_m[pm+1]-proc_start_m[pm];
        auto bs    = proc_start_n[pn+1]-proc_start_n[pn];
        recv_layouts[iproc] = proc_inner_layout(iproc,global_n2b());
        recv_displs[iproc] = proc_start_m[pm]*n_global + proc_start_n[pm];
      }
    //   // comm.gatherv
    //   //   ( 0,
    //   //           first_data_point,send_layout,
    //   //     gathered_data->data(),recv_layouts,recv_displs );
    //   // view as zero-border array
    //   bordered_array_seq<real>( m_global,n_global,gathered_data->data() ).view(caption);
    // } else {
    //   comm.reduce( mpl::plus<size_t>(),0,local_size );
      comm.gatherv( 0, first_data_point,send_layout );
    }
  };

  /*
   * Logging
   */
  template< typename real >
  void distributed_array<real>::log_flops( float n ) {
    subdomain->log_flops(n); };
  template< typename real >
  void distributed_array<real>::log_bytes( float n ) {
    subdomain->log_bytes(n); };
  template< typename real >
  std::pair<float,float> distributed_array<real>::log_report() {
    // TODO this needs a reduction
    return subdomain->log_report(); };

    // upcast base to derived type
    // const auto& other = dynamic_cast<const linalg::distributed_array<real>&>(_other);

};

namespace linalg {
  template class distributed_array<float>;
  template class distributed_array<double>;
};

