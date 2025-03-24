/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2025 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** main.cpp : general main for diff2d codes
 ****
 ****************************************************************/

X->set_value( 1.,trace );
if (view)
  X->view("Original");

auto xnorm = X->l2norm();
if ( trace and procno==0 )
  std::cout << std::format("x norm: {}\n",xnorm);

using myclock = std::chrono::steady_clock;
auto start_time = myclock::now();

//codesnippet d2dmain
for ( int it=0; it<itcount; it++ ) {
  /*
   * Matrix-vector product
   */
  Y->central_difference_from( *X,trace );
  //omitbegin
  if (view)
    Y->view("Operator applied");
  //omitend
  // norm computation
  auto bnorm = Y->l2norm();
  //omitbegin
  if ( trace and procno==0 )
    std::cout << std::format("[{:>2}] y norm: {}\n",it,bnorm);
  //omitend
  // scale
  X->scale_interior( *Y,1./bnorm );
  //omitbegin
  if ( view )
    X->view("Scaled");
  //omitend
 }
//codesnippet end

/*
 * Time reporting
 */
auto duration = myclock::now()-start_time;
auto millisec_duration = 
  std::chrono::duration_cast<std::chrono::microseconds>(duration)/1000;
auto msec = millisec_duration.count();
if ( procno==0 )
  std::cout << std::format("Time: {:>6} msec\n",msec);

/*
 * Other stats
 */
auto [flops,bytes] = X->log_report();
if ( procno==0 )
  std::cout << std::format("{}-BW: (estimated Gb/sec) {:7.3}\n",
        prefix,bytes/(msec/1.e3)*1.e-9);
if ( procno==0 )
  std::cout << std::format("{}-Flops: {:7.3}\n",
        prefix,flops/(msec/1.e3)*1.e-9);
