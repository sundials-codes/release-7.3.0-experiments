cmake_minimum_required(VERSION 3.22)
project(gray_scott CXX C)

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
  GIT_TAG v7.2.1
  GIT_SHALLOW TRUE
)
set(ENABLE_OPENMP ON)
set(BUILD_CVODE OFF)
set(BUILD_CVODES OFF)
set(BUILD_IDA OFF)
set(BUILD_IDAS OFF)
set(BUILD_KINSOL OFF)
set(SUNDIALS_LOGGING_LEVEL 0)
set(SUNDIALS_INDEX_SIZE 32)
set(EXAMPLES_ENABLE_C OFF)
set(BUILD_SHARED_LIBS OFF)
FetchContent_MakeAvailable(sundials)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC
  SUNDIALS::core
  SUNDIALS::arkode
  SUNDIALS::nvecopenmp
  OpenMP::OpenMP_CXX)
