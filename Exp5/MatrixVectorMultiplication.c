#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int my_rank, comm_sz;

void ReadMatrix(long m, long n, double* local_mat_cols){
    double* matrix;
    long i;

    if(my_rank == 0){
        matrix = (double*)malloc(m * n * sizeof(double));
        for(i = 0; i < m * n; i++){
            scanf("%lf", matrix + i);
        }

    free(matrix);
    }else{

    }
}

void ReadVector(){

}

int main(int argc, char** argv){
    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
}
