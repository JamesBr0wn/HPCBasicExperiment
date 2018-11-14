cmake_minimum_required(VERSION 3.10)
project(HPCBasicExperiment)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(MPI_Proc_Num 4)

add_subdirectory(Exp1 Exp1)
add_subdirectory(Exp2 Exp2)
add_subdirectory(Exp3 Exp3)
add_subdirectory(Exp4 Exp4)