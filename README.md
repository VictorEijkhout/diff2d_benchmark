# Diff2d benchmark

## compilation

Go into `code/diff2d` and use make or cmake:

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

