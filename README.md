# Lattice reduction - Benchmarking LLL vs Jacobi Reduction in various settings

## Compiler and Environment



## Libraries 

- [`The new Number Theory Library`](http://www.prism.uvsq.fr/~gama/newntl.html) is used for the Jacobi method implementation. Note: newNTL needs to be built with GCC and does not build with Clang. When G++ is installed (`brew install gcc`): `make CXX=/usr/local/bin/g++-(version)`
- [`fplll`](http://perso.ens-lyon.fr/damien.stehle/fplll/) is used as a reference LLL implementation for benchmarking. 