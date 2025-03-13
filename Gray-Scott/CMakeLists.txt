cmake_minimum_required(VERSION 3.22)
project(gray_scott CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-Wall -Wextra -Wpedantic -Werror)

find_package(OpenMP REQUIRED)
# TODO: automatically download specific version of SUNDIALS
find_package(SUNDIALS REQUIRED CONFIG PATHS "../../sundials/install/lib64/cmake/sundials/")
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC
  SUNDIALS::core
  SUNDIALS::arkode
  SUNDIALS::nvecopenmp
  OpenMP::OpenMP_CXX)
