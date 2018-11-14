#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv){
    int my_id, comm_sz;
    long total_hits, hits_in_cycle, my_total_hits, my_hits_in_cycle, i;
    double pi, x, y, dist_sq;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    srand((unsigned int)(time(NULL)) + my_id);

    if(my_id == 0){
        // 进程0向其他进程广播总投掷次数
        total_hits = 2000000;
        for(i = 1; i < comm_sz; i++){
            MPI_Send(&total_hits, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
        }
    }else{
        // 其他进程接收总投掷次数
        MPI_Status status;
        MPI_Recv(&total_hits, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
    }

    // 蒙特卡洛方法估计Pi值
    my_total_hits = total_hits / comm_sz;
    my_hits_in_cycle = 0;
    for(i = 0; i < my_total_hits; i++){
        x = rand() / (double)RAND_MAX;
        y = rand() / (double)RAND_MAX;
        dist_sq = x * x + y * y;
        if(dist_sq <= 1){
            my_hits_in_cycle++;
        }
    }
    double my_pi = 4 * my_hits_in_cycle / (double)my_total_hits;
    printf("Pi of process %d:\t%.4f\n", my_id, my_pi);

    // MPI_Reduce合并各进程结果
    MPI_Reduce(&my_hits_in_cycle, &hits_in_cycle, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if(my_id == 0){
        pi = 4 * hits_in_cycle / (double)total_hits;
        printf("Pi Estimate:\t%.4f\n", pi);
    }
    MPI_Finalize();
    return 0;
}