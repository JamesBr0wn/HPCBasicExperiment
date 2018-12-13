cmake_minimum_required(VERSION 3.10)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

list(APPEND ExecutableList CountSort)

foreach(executable IN LISTS ExecutableList)
    add_executable(${executable} ${executable}.c)
endforeach()