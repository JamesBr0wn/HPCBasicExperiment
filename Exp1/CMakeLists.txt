cmake_minimum_required(VERSION 3.10)

list(APPEND ExecutableList NoRotate)
list(APPEND ExecutableList RotateNaive)
list(APPEND ExecutableList RotateSummary)
list(APPEND ExecutableList RotateVerOne)
list(APPEND ExecutableList RotateVerTwo)
list(APPEND ExecutableList RotateVerThree)

foreach(executable IN LISTS ExecutableList)
    add_executable(${executable} ${executable}.cpp)
endforeach()