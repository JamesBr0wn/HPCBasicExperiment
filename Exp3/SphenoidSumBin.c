#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv){
    // 变量声明
    int my_id, comm_sz, data, temp_data, comm_count = 1, dst;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // srand((unsigned int)(time(NULL)) + my_id);
    // data = min + (max - min) * rand() / (double)RAND_MAX;
    data = my_id + 1;

    while(comm_count * 2 <= comm_sz){
        comm_count *= 2;    // 计数值乘2，进入下一轮通信
        if(my_id % comm_count < comm_count / 2){
            dst = my_id + comm_count / 2;
        }else{
            dst = my_id - comm_count / 2;
        }
        // 配对进程交换部分和
        MPI_Sendrecv(&data, 1, MPI_INT, dst, comm_count, &temp_data, 1, MPI_INT, dst, comm_count, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d exchange local data with %d\n", my_id, dst);
        data += temp_data;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    printf("Total sum of prcess %d is:\t%d\n", my_id, data);
    MPI_Finalize();
    return 0;
}