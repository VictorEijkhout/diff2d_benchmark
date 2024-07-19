# Diff2d benchmark

This is a benchmark for parallel computation of a "power method":
stencil-based matrix-vector product, norm, scaling, and repeat that.

So far there are C++ implementations based on OpenMP, Kokkos, Sycl, MPI.
Contributed implementations welcome.

## Compilation

Go into `code/diff2d`. Calling `make` without arguments
tells you all the make rules. For compilation use make or cmake:

```
make compile VARIANTS="seq oned"
```

or

```
make cmake VARIANTS="kokkos sycl"
```

You can of course run cmake outside of make:

```
variant=span
mkdir build
ln cmake/CMakeLists.txt $variant
cmake -B build -S $variant -D VARIANT=$variant
```

Until I figure out what targets are provided in `mdspan` from Kokkos,
do additionally

```
cmake \
    -D MDSPAN_INC=${TACC_MDSPAN_INC}
```

