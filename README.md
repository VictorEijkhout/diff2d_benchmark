# Diff2d benchmark

This is a benchmark for parallel computation of a "power method":
stencil-based matrix-vector product, norm, scaling, and repeat that.

So far there are C++ implementations based on OpenMP, Kokkos, Sycl, MPI.
Contributed implementations welcome.

## Compilation

### Prerequisites

This software uses the package [cxxopts](https://github.com/jarro2783/cxxopts)
and [mdspan](https://github.com/kokkos/mdspan).

If you want to install those yourself:

 - add the `.pc` files from `cxxopts` to the `PKG_CONFIG_PATH`
 - add the `mdspan` installation directory to the `CMAKE_PREFIX_PATH`.

### Makefile compilation

Go into `code/diff2d`. Calling `make` without arguments
tells you all the make rules. For compilation use make or cmake:

```
make compile VARIANTS="seq oned"
```

Set the variable `TACC_MDSPAN_INC`
to the location of the header files.

### CMake compilation

Go into `code/diff2d`. Drive the cmake installation with the makefile:

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

using of course your own actual installation path.
