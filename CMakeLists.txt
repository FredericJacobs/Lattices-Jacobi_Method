cmake_minimum_required(VERSION 2.8.4)
project(JacobiMethod)

SET (CMAKE_CXX_COMPILER "/usr/local/bin/g++-4.9")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include /usr/local/include/newNTL)

find_library(NEWNTL newntl)
find_library(MPFR   mpfr)
find_library(GMP    gmp)
find_library(FPLLL  fplll)

set(SOURCE_FILES main.cpp src/JacobiMethod.cpp)
add_executable(JacobiMethod ${SOURCE_FILES})
target_link_libraries(JacobiMethod ${NEWNTL} ${MPFR} ${GMP} ${FPLLL})


