cmake_minimum_required(VERSION 3.31)
project(Diploma_computing)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# путь к установленному libomp
set(LIBOMP_ROOT "/opt/homebrew/opt/libomp")
set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_C_LIB_NAMES "omp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "${LIBOMP_ROOT}/lib/libomp.dylib")
set(OpenMP_omp_LIBRARY_DEBUG "${OpenMP_omp_LIBRARY}")

find_package(OpenMP REQUIRED)


add_executable(Diploma_computing
        crochemore_paralell.cpp
)


target_include_directories(Diploma_computing PRIVATE ${LIBOMP_ROOT}/include)
target_link_directories(Diploma_computing PRIVATE ${LIBOMP_ROOT}/lib)
target_link_libraries(Diploma_computing PRIVATE OpenMP::OpenMP_CXX)