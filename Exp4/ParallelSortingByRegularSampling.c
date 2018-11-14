#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <mpi.h>

int main(int argc, char** argv){
    // 变量声明
    int my_rank, comm_sz;
    size_t init_data_count, data_count, local_data_count;
    long i, *data, *local_data;

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
        fread(&data, sizeof(long), init_data_count - 1, input_file);

        for(i = data_count - 1; i >= init_data_count; i--){
            data[i] = LONG_MAX;
        }
    }
    return 0;
}
