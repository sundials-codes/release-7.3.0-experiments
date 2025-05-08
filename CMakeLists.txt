cmake_minimum_required(VERSION 3.22)
project(sundials_release_730_paper CXX C)

set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_options(-fno-math-errno -march=native)

find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
  sundials
  GIT_REPOSITORY https://github.com/LLNL/sundials
  GIT_TAG v7.3.0
  GIT_SHALLOW TRUE
)
set(ENABLE_OPENMP ON)
set(BUILD_CVODE OFF)
set(BUILD_CVODES ON)
set(BUILD_IDA OFF)
set(BUILD_IDAS OFF)
set(BUILD_KINSOL OFF)
set(SUNDIALS_LOGGING_LEVEL 0)
set(SUNDIALS_INDEX_SIZE 32)
set(EXAMPLES_ENABLE_C OFF)
set(BUILD_SHARED_LIBS OFF)
FetchContent_MakeAvailable(sundials)

add_executable(gray_scott Gray-Scott/main.cpp)
target_link_libraries(gray_scott PUBLIC
  SUNDIALS::core
  SUNDIALS::arkode
  SUNDIALS::nvecopenmp
  OpenMP::OpenMP_CXX)

add_executable(lotka_volterra_arkode Lotka-Volterra/lotka_volterra_arkode.c)
target_link_libraries(lotka_volterra_arkode PUBLIC
  SUNDIALS::core
  SUNDIALS::arkode)

add_executable(lotka_volterra_cvodes Lotka-Volterra/lotka_volterra_cvodes.c)
  target_link_libraries(lotka_volterra_cvodes PUBLIC
    SUNDIALS::core
    SUNDIALS::cvodes)
