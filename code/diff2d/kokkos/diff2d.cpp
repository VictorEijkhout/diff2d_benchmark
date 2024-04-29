#include <iostream>
#include <chrono>
#include <Kokkos_Core.hpp>

#define ENABLE_CUDA false

#include "../lib/options.hpp"
using real = float;

int main(int argc, char *argv[]) {

  auto [exit,msize,nsize,border,itcount,trace,view] =
    parse_options(argc,argv,"Kokkos version using 2D loop");
  if (exit) return 0;
  const std::string prefix{"kks"};
  const int procno{0};

  Kokkos::initialize(argc, argv);
  {
    using MemSpace = Kokkos::HostSpace;  // Simplify to HostSpace for printing
    using Layout = Kokkos::LayoutRight;

    using ViewMatrixType = Kokkos::View<real**, Layout, MemSpace>;
    ViewMatrixType x("x", msize,nsize ), Ax("Ax", msize,nsize);

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

    for (int repeat = 0; repeat<itcount; ++repeat) {
      Kokkos::parallel_for
	("StencilApplication",
	 Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {msize-1,nsize-1}),
	 KOKKOS_LAMBDA(int i, int j) {
	  Ax(i, j) = 4 * x(i, j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1);
	});

      //Kokkos::deep_copy(h_Ax, Ax);

      // Compute and apply scaling
      real norm = 0.0;
      Kokkos::parallel_reduce
	("Compute norm",
	 Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {msize-1, nsize-1}),
	 KOKKOS_LAMBDA(int i, int j, real& update) {
	  update += Ax(i, j) * Ax(i, j);
	}, norm);
      norm = std::sqrt(norm);

      Kokkos::parallel_for
	("Update x",
	 Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {msize-1, nsize-1}),
	 KOKKOS_LAMBDA(int i, int j) {
	  x(i, j) = Ax(i, j) / norm;
	});
      
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
