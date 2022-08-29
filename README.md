# Hot plasma kernel

Hot plasma is characterised by a non-local conductivity kernel.
This code implements a simple finite element model, in 1D, in an annulus geometry.
The required integrals are evaluated, and a routine is provided to rapidly compute matrix entries.

# Installation
You need a C++ compiler, `cmake` and `make`.

Then compile:

` mkdir build ; cd build`

`cmake .. -DCMAKE_BUILD_TYPE=Release`

`make -j NTHREADS`

(e.g. `make -j 6` to compile with 6 threads)

# License
HierarchicalRF  Copyright (C) 2022  Mike Machielsen\
This program comes with ABSOLUTELY NO WARRANTY.\
This is free software, and you are welcome to redistribute it\
under certain conditions; see license for details.