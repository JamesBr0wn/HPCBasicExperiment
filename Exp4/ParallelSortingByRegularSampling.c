#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <mpi.h>

int cmp ( const void *a , const void *b ) {
    return *(long*)a > *(long*)b ? 1 : -1;
}

int main(int argc, char** argv){
    // 变量声明
    int my_rank, comm_sz;
    size_t init_data_count, data_count, local_data_count;
    long i, *data = NULL, *local_data = NULL, *sample = NULL, *local_sample = NULL;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // 数据分发
    if(my_rank == 0){
        FILE* input_file = fopen("psrs_data", "rb");
        fread(&init_data_count, sizeof(long), 1, input_file);
        data_count = init_data_count % comm_sz == 0 ? init_data_count : init_data_count  + comm_sz - init_data_count % comm_sz;

        data = (long*)malloc(data_count * sizeof(long));
        fread(data, sizeof(long), init_data_count, input_file);
        for(i = data_count - 1; i >= init_data_count; i--){
            data[i] = LONG_MAX;
        }

        fclose(input_file);

        local_data_count = data_count / comm_sz;
        local_data = (long*)malloc(local_data_count * sizeof(long));

        MPI_Bcast(&local_data_count, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Scatter(data, (int)local_data_count, MPI_LONG, local_data, (int)local_data_count, MPI_LONG, 0, MPI_COMM_WORLD);
    }else{
        MPI_Bcast(&local_data_count, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        local_data = (long*)malloc(local_data_count * sizeof(long));
        MPI_Scatter(data, (int)local_data_count, MPI_LONG, local_data, (int)local_data_count, MPI_LONG, 0, MPI_COMM_WORLD);
    }

    // 本地排序
    qsort(local_data, local_data_count, sizeof(local_data[0]), cmp);

    printf("*%d %ld\n", my_rank, local_data_count);
    for(i = 0; i < local_data_count; i++){
        printf("-%d %ld \n", my_rank, local_data[i]);
    }
    putchar('\n');

    MPI_Finalize();
    return 0;
}
