# Hot plasma kernel

Hot plasma is characterised by a non-local conductivity kernel.
This code implements a simple finite element model, in 1D, in an annulus geometry.
The required integrals are evaluated, and a routine is provided to rapidly compute matrix entries.

# Installation
You need a C++ compiler, `cmake`, `make` and the Boost (version>=1.70 because of the complex quadrature) and Eigen3 libraries.
You need to define an environment variable `BOOST_INCLUDE_DIR` equal to the path in which you extracted the Boost.
Similarly, for cmake to be able to find Eigen, define `EIGEN3_INCLUDE_DIR` to be equal to the path in which eigen is installed.
E.g. `export EIGEN3_INCLUDE_DIR=$HOME/Downloads/eigen/`. You also need OpenMP, but most compilers these days come with this included.
Lastly, you also need the GNU scientific library.

Then compile:

` mkdir build ; cd build`

`cmake .. -DCMAKE_BUILD_TYPE=Release`

`make -j NTHREADS`

(e.g. `make -j 6` to compile with 6 threads)

# Running the program
The executable in inside the build directory.

`./hotPlasmaKernel.x`

You can control the number of threads by setting the environment variable `OMP_NUM_THREADS`.

# Misc

Optionally, you can generate a documentation file with doxygen by running `doxygen` in the root dir (i.e. where the Doxyfile is located).
Then open html/index.html to read it. To generate the pdf version, navigate to `latex/` and use `make`, this requires `pdflatex`.

# License
HierarchicalRF  Copyright (C) 2022  Mike Machielsen\
This program comes with ABSOLUTELY NO WARRANTY.\
This is free software, and you are welcome to redistribute it\
under certain conditions; see license for details.