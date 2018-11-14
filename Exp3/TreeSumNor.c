#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv){
    // 变量声明
    int my_id, comm_sz, data, temp_data, comm_count = 1;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // srand((unsigned int)(time(NULL)) + my_id);
    // data = min + (max - min) * rand() / (double)RAND_MAX;
    data = my_id + 1;

    // 聚集
    int max_count = pow(2, (int)ceil(log(comm_sz) * 1.0 / log(2)));
    MPI_Status status;
    while(comm_count < max_count){
        comm_count *= 2;    // 计数值乘2，进入下一轮通信
        if(my_id % comm_count != 0){    // 对于发送进程，将自己对应的部分和发送出去并跳出
            printf("Process %d send local data to %d\n", my_id, my_id / comm_count * comm_count);
            MPI_Send(&data, 1, MPI_INT, my_id / comm_count * comm_count, comm_count, MPI_COMM_WORLD);
            break;
        }else{  // 对于接收进程，接收来自其他进程的部分和，与自己的部分和累加
            if(my_id + comm_count / 2 < comm_sz){
                printf("Process %d receieve local data from %d\n", my_id, my_id + comm_count / 2);
                MPI_Recv(&temp_data, 1, MPI_INT, my_id + comm_count / 2, comm_count, MPI_COMM_WORLD, &status);
                data += temp_data;
            }
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // 分发
    while(comm_count / 2 > 0){
        comm_count /= 2;
        if(my_id % (comm_count * 2) == 0){
            if(my_id + comm_count < comm_sz){   // 对于发送进程，将全局和发送出去
                printf("Process %d send global data to %d\n", my_id, my_id + comm_count);
                MPI_Send(&data, 1, MPI_INT, my_id + comm_count, comm_count, MPI_COMM_WORLD);
            }
        }else if(my_id % comm_count == 0){  // 对于接收进程，接收来自其他进程的全局和
            MPI_Recv(&data, 1, MPI_INT, MPI_ANY_SOURCE, comm_count, MPI_COMM_WORLD, &status);
            printf("Process %d receieve global data from %d\n", my_id, status.MPI_SOURCE);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("Total sum of prcess %d is:\t%d\n", my_id, data);
    MPI_Finalize();
    return 0;
}