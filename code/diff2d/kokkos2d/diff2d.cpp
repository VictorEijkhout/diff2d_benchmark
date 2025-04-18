#include <iostream>
#include <chrono>
#include <Kokkos_Core.hpp>

#define ENABLE_CUDA false

#include "../lib/options.hpp"
using real = float;

int main(int argc, char *argv[]) {

  auto [exit,msize,nsize,border,itcount,gpu,trace,view] =
    parse_options(argc,argv,"Kokkos version using 2D loop");
  if (exit) return 0;
  const std::string prefix{"kks"};
  const int procno{0};

  Kokkos::initialize(argc, argv);
  {
    //codesnippet kokkosbufcreate
    using MemSpace = Kokkos::HostSpace;
    using Layout = Kokkos::LayoutRight;
    using HostMatrixType = Kokkos::View<real**, Layout, MemSpace>;
    HostMatrixType x("x", msize,nsize );
    //codesnippet end
    HostMatrixType Ax("Ax", msize,nsize);

    // Initialize matrix with boundaries
#if 0
    Kokkos::parallel_for
      ("Initialize matrix",
       Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {msize,nsize}),
       KOKKOS_LAMBDA(int i, int j) {
        if (i > 0 && i < msize-1 && j > 0 && j < nsize-1) {
          x(i, j) = 1.0;  // Interior
        } else {
          x(i, j) = 0.0;  // Boundary
        }
      });
#else
    for (int i=0; i<msize; i++)
      for (int j=0; j<nsize; j++)
        if (i > 0 && i < msize-1 && j > 0 && j < nsize-1) {
          x(i, j) = 1.0;  // Interior
        } else {
          x(i, j) = 0.0;  // Boundary
        }
#endif

    Kokkos::View<real**, Kokkos::HostSpace> h_x = Kokkos::create_mirror_view(x);
    Kokkos::View<real**, Kokkos::HostSpace> h_Ax = Kokkos::create_mirror_view(Ax);
    Kokkos::deep_copy(h_x, x);
    //        printMatrix("Original matrix:", h_x, N);

    using myclock = std::chrono::steady_clock;
    auto start_time = myclock::now();

    for (int it = 0; it<itcount; ++it) {
      Kokkos::parallel_for
        ("StencilApplication",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {msize-1,nsize-1}),
         KOKKOS_LAMBDA(int i, int j) {
          Ax(i, j) = 4 * x(i, j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1);
        });

      //Kokkos::deep_copy(h_Ax, Ax);

      // Compute and apply scaling
      //codesnippet kokkosreduce
      real norm = 0.0;
      Kokkos::parallel_reduce
        ("Compute norm squared",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {msize-1, nsize-1}),
         KOKKOS_LAMBDA(int i, int j, real& accum) {
          accum += Ax(i, j) * Ax(i, j);
        }, norm);
      norm = std::sqrt(norm);
      //codesnippet end

      if ( trace and procno==0 )
        std::cout << std::format("[{:>2}] y norm: {}\n",it,norm);

      //codesnippet kokkosbufaccess
      Kokkos::parallel_for
        ("Update x",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>
             ({1, 1}, {msize-1, nsize-1}),
         KOKKOS_LAMBDA(int i, int j) {
          x(i, j) = Ax(i, j) / norm;
        });
      //codesnippet end
      
      //Kokkos::deep_copy(h_x, x);
    }

    /*
     * Time reporting
     */
    auto duration = myclock::now()-start_time;
    auto millisec_duration = 
      std::chrono::duration_cast<std::chrono::microseconds>(duration)/1000;
    auto msec = millisec_duration.count();
    if ( procno==0 )
      std::cout << std::format("Time: {:>6} msec\n",msec);
  }

  Kokkos::finalize();
  return 0;
}
