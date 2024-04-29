#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <Kokkos_Core.hpp>

#define ENABLE_CUDA false

void checkArgs(int &n, int &nrepeat);
/*
void printMatrix(const char* title, Kokkos::View<double**, Kokkos::HostSpace> matrix, int N) {
    std::cout << title << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
*/
int main(int argc, char *argv[]) {
    int n = 4;  // Internal grid size
    int N = n + 2;  // Total size including boundaries
    int nrepeat = 100;  // Reduce to 1 for demonstration purposes

    checkArgs(n, nrepeat);

    Kokkos::initialize(argc, argv);
    {
        using MemSpace = Kokkos::HostSpace;  // Simplify to HostSpace for printing
        using Layout = Kokkos::LayoutRight;

        using ViewMatrixType = Kokkos::View<double**, Layout, MemSpace>;
        ViewMatrixType x("x", N, N), Ax("Ax", N, N);

        // Initialize matrix with boundaries
        Kokkos::parallel_for("Initialize matrix", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {N, N}),
                             KOKKOS_LAMBDA(int i, int j) {
                                 if (i > 0 && i < N-1 && j > 0 && j < N-1) {
                                     x(i, j) = 1.0;  // Interior
                                 } else {
                                     x(i, j) = 0.0;  // Boundary
                                 }
                             });

        Kokkos::View<double**, Kokkos::HostSpace> h_x = Kokkos::create_mirror_view(x);
        Kokkos::View<double**, Kokkos::HostSpace> h_Ax = Kokkos::create_mirror_view(Ax);
        Kokkos::deep_copy(h_x, x);
//        printMatrix("Original matrix:", h_x, N);

        for (int repeat = 0; repeat < nrepeat; ++repeat) {
            Kokkos::parallel_for("StencilApplication", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {N-1, N-1}),
                                 KOKKOS_LAMBDA(int i, int j) {
                                     Ax(i, j) = 4 * x(i, j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1);
                                 });

            Kokkos::deep_copy(h_Ax, Ax);
//            printMatrix("After operator applied:", h_Ax, N);

            // Compute and apply scaling
            double norm = 0.0;
            Kokkos::parallel_reduce("Compute norm", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {N-1, N-1}),
                                    KOKKOS_LAMBDA(int i, int j, double& update) {
                                        update += Ax(i, j) * Ax(i, j);
                                    }, norm);
            norm = std::sqrt(norm);

            Kokkos::parallel_for("Update x", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {N-1, N-1}),
                                 KOKKOS_LAMBDA(int i, int j) {
                                     x(i, j) = Ax(i, j) / norm;
                                 });

            Kokkos::deep_copy(h_x, x);
//            printMatrix("After scaling:", h_x, N);
        }

    }

        Kokkos::finalize();
    return 0;
}

void checkArgs(int &n, int &nrepeat) {
    if (n < 3 || nrepeat < 1) {
        std::cerr << "Error: grid size n must be at least 3 and nrepeat must be greater than 0." << std::endl;
        exit(1);
    }
}

