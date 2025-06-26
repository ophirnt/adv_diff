cmake_minimum_required(VERSION 3.10)
project(main CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# PETSc location
set(PETSC_DIR "$ENV{HOME}/petsc")
set(PETSC_ARCH "arch-linux-c-opt")
set(PETSC_PATH "${PETSC_DIR}/${PETSC_ARCH}")

# Include and library paths
include_directories(${PETSC_PATH}/include ${PETSC_DIR}/include)
link_directories(${PETSC_PATH}/lib)

add_executable(main src/main.cpp)

# Link against PETSc and common system dependencies
target_link_libraries(main
    petsc
    m       # math
    pthread
    X11     # if needed
)
