cmake_minimum_required(VERSION 3.10)
project(main CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find PETSc
find_package(PETSc REQUIRED)

include_directories(${PETSC_INCLUDES} src)

add_executable(main src/main.cpp)

target_link_libraries(main ${PETSC_LIBRARIES}) 