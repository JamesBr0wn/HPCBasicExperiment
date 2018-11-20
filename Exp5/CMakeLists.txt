cmake_minimum_required(VERSION 3.10)

list(APPEND ExecutableList MatrixVectorMultiplication)

foreach(executable IN LISTS ExecutableList)
    add_executable(${executable} ${executable}.c)
    target_link_libraries(${executable} ${MPI_C_LIBRARIES} m)

    if(MPI_COMPILE_FLAGS)
        set_target_properties(${executable} PROPERTIES
                COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(${executable} PROPERTIES
                LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()

    #    add_custom_command(TARGET ${executable} POST_BUILD COMMAND time mpirun -np ${MPI_Proc_Num} ./${executable} > ${executable}.log 2>> ${executable}.log)
endforeach()