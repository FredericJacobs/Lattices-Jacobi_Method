cmake_minimum_required(VERSION 2.8.4)


project(JacobiMethodTests)

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

set(SOURCE_FILES src/JacobiTests.cpp src/JacobiMethod.cpp src/MatrixFactory.cpp src/ReductionQualityChecker.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_executable(JacobiMethodTests ${SOURCE_FILES})
target_link_libraries(JacobiMethodTests ${NEWNTL} ${MPFR} ${GMP} ${FPLLL})


project(LLLThenJacobi)

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

set(SOURCE_FILES src/LLLThenJacobi.cpp src/JacobiMethod.cpp src/MatrixFactory.cpp src/ReductionQualityChecker.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_executable(LLLThenJacobi ${SOURCE_FILES})
target_link_libraries(LLLThenJacobi ${NEWNTL} ${MPFR} ${GMP} ${FPLLL})
