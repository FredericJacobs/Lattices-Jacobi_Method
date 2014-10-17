# Building

This projects uses [CMAKE](http://www.cmake.org/) as a cross-platform build processing. The code of the project is portable, but needs to be linked against newNTL which only builds with GCC.

Make sure to export the CC and CXX GCC environment variables before running CMAKE.

On a Mac, that has a Homebrew installed (`brew install gcc`) version of GCC, that can be done by running the following export commands:

```bash
export CC=/usr/local/bin/gcc-4.9
export CXX=/usr/local/bin/g++-4.9
```

Running the `cmake CMakeLists.txt` command should generate a Makefile for your system.

Simply run `make` after that to build the targets.