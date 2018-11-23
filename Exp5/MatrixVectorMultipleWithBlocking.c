#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int my_rank, my_row_rank, my_col_rank, my_pivot_rank, comm_sz, row_comm_sz, col_comm_sz, pivot_comm_sz;
MPI_Comm total_comm, row_comm, col_comm, pivot_comm;

void create_comm(){
    int sub_comm_sz = (int)sqrt(comm_sz);
    int ndims = 2;
    int dims[2] = {sub_comm_sz, sub_comm_sz};
    int periods[2] = {0, 0};
    int reorder = 0;
    int remain_dims[2] = {1, 1};
    int pivots[sub_comm_sz];
    MPI_Group total_group, pivot_group;
    int i;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &total_comm);

    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(total_comm, remain_dims, &row_comm);
    MPI_Comm_size(row_comm, &row_comm_sz);
    MPI_Comm_rank(row_comm, &my_row_rank);

    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(total_comm, remain_dims, &col_comm);
    MPI_Comm_size(col_comm, &col_comm_sz);
    MPI_Comm_rank(col_comm, &my_col_rank);

    for(i = 0; i < sub_comm_sz; i++){
        pivots[i] = i * (sub_comm_sz + 1);
    }

    MPI_Comm_group(MPI_COMM_WORLD, &total_group);
    MPI_Group_incl(total_group, sub_comm_sz, pivots, &pivot_group);
    MPI_Comm_create(MPI_COMM_WORLD, pivot_group, &pivot_comm);
    if(my_row_rank == my_col_rank){
        MPI_Comm_size(pivot_comm, &pivot_comm_sz);
        MPI_Comm_rank(pivot_comm, &my_pivot_rank);
    }else{
        pivot_comm_sz = sub_comm_sz;
        my_pivot_rank = -1;
    }
}

void get_size(int* m_ptr, int* local_m_ptr, int* n_ptr, int* local_n_ptr){
    if(my_rank == 0){
        printf("Please enter the size of matrix and vector:\t");
        scanf("%d", n_ptr);
    }
    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    *m_ptr = *n_ptr;
    *local_m_ptr = *m_ptr / row_comm_sz;
    *local_n_ptr = *n_ptr / col_comm_sz;
}

int main(int argc, char** argv){
    int m, n, local_m, local_n;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    create_comm();

    get_size(&m, &local_m, &n, &local_n);

    printf("%2d %2d %2d %2d\n", my_rank, my_row_rank, my_col_rank, my_pivot_rank);

    MPI_Finalize();

    return 0;
}