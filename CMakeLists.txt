cmake_minimum_required(VERSION 2.8.4)
project(JacobiMethod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include /usr/local/include/newNTL)

find_library(NEWNTL newntl)
find_library(MPFR   mpfr)
find_library(GMP    gmp)
find_library(FPLLL  fplll)

set(SOURCE_FILES main.cpp src/JacobiMethod.cpp src/MatrixFactory.cpp)
add_executable(JacobiMethod ${SOURCE_FILES})
target_link_libraries(JacobiMethod ${NEWNTL} ${MPFR} ${GMP} ${FPLLL})
