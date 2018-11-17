#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <mpi.h>

int cmp ( const void *a , const void *b ) {
    return *(long*)a > *(long*)b ? 1 : -1;
}

void read_from_file(const char* file_name, int comm_sz, size_t* init_data_count, size_t* data_count, long** data){
    size_t i;
    FILE* input_file = fopen(file_name, "rb");
    fread(init_data_count, sizeof(long), 1, input_file);

    *data_count = (*init_data_count) % comm_sz == 0 ? (*init_data_count) : (*init_data_count)  + comm_sz - (*init_data_count) % comm_sz;

    *data = (long*)malloc((*data_count) * sizeof(long));
    fread(*data, sizeof(long), *init_data_count, input_file);

    for(i = *data_count - 1; i >= *init_data_count; i--){
        (*data)[i] = LONG_MAX;
    }
    fclose(input_file);
}

void deliver_data(const size_t data_count, const long *data, int my_rank, int comm_sz, size_t *local_data_count, long **local_data){
    if(my_rank == 0) {
        *local_data_count = data_count / comm_sz;
    }
    MPI_Bcast(local_data_count, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    *local_data = (long*)malloc((*local_data_count) * sizeof(long));
    MPI_Scatter(data, (int)(*local_data_count), MPI_LONG, *local_data, (int)(*local_data_count), MPI_LONG, 0, MPI_COMM_WORLD);
}

void select_pivot(size_t local_data_count, const long* local_data, int my_rank, int comm_sz, long** pivots){
    long *sample = NULL, *local_sample = NULL;
    size_t step = (local_data_count + comm_sz - 1)/ comm_sz, i = 0;
    local_sample = (long*)malloc(comm_sz * sizeof(long));
    for(i = 0; i < comm_sz; i++){
        local_sample[i] = local_data[i * step];
    }

    sample = (long*)malloc(comm_sz * comm_sz * sizeof(long));
    MPI_Gather(local_sample, comm_sz, MPI_LONG, sample, comm_sz, MPI_LONG, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
        qsort(sample, (size_t)comm_sz * comm_sz, sizeof(long), cmp);
    }
    free(local_sample);

    *pivots = (long*)malloc((comm_sz - 1) * sizeof(long));
    if(my_rank == 0){
        for(i = 0; i < comm_sz - 1; i++){
            (*pivots)[i] = sample[(i+1) * comm_sz + comm_sz / 2 -1];
        }
    }
    MPI_Bcast(*pivots, comm_sz - 1, MPI_LONG, 0, MPI_COMM_WORLD);
    free(sample);
}

void exchange_data(size_t local_data_count, long *local_data, long *pivots, int comm_sz, size_t *exchanged_data_count, long **exchanged_data, int **merge_counts){
    int i = 0, j = 0, *send_counts, *recv_counts, *send_disps, *recv_disps;

    send_counts = (int*)malloc(comm_sz * sizeof(int));
    recv_counts = (int*)malloc(comm_sz * sizeof(int));

    for(i = 0; i < comm_sz; i++){
        send_counts[i] = 0;
        recv_counts[i] = 0;
    }

    for(i = 0; i < local_data_count; i++){
        while(j < comm_sz - 1 && local_data[i] > pivots[j]){
            j++;
        }
        if(j == comm_sz - 1){
            send_counts[comm_sz - 1] = (int)local_data_count - i;
            break;
        }
        send_counts[j]++;
    }

    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

    for(i = 0; i < comm_sz; i++){
        *exchanged_data_count += recv_counts[i];
    }
    *exchanged_data = (long*)malloc((*exchanged_data_count) * sizeof(long));

    send_disps = (int*)malloc(comm_sz * sizeof(int));
    recv_disps = (int*)malloc(comm_sz * sizeof(int));

    send_disps[0] = 0;
    recv_disps[0] = 0;
    for(i = 1; i < comm_sz; i++){
        send_disps[i] = send_disps[i-1] + send_counts[i-1];
        recv_disps[i] = recv_disps[i-1] + recv_counts[i-1];
    }

    MPI_Alltoallv(local_data, send_counts, send_disps, MPI_LONG, *exchanged_data, recv_counts, recv_disps, MPI_LONG, MPI_COMM_WORLD);
    *merge_counts = recv_counts;

    free(pivots);
    free(local_data);
    free(send_counts);
    free(send_disps);
    free(recv_disps);
}

void merge_data(size_t exchanged_data_count, long* exchanged_data, int comm_sz, int* merge_counts, long** merged_data){
    int i = 0, j = 0, temp = 0, *merge_indexs, *merge_ends;
    merge_indexs = (int*)malloc(comm_sz * sizeof(int));
    merge_ends = (int*)malloc(comm_sz * sizeof(int));
    for(i = 0; i < comm_sz; i++){
        merge_indexs[i] = temp;
        temp += merge_counts[i];
        merge_ends[i] = temp;
    }

    *merged_data = (long*)malloc(exchanged_data_count * sizeof(long));
    long min;
    int min_idx;
    for(i = 0; i < exchanged_data_count; i++){
        min = LONG_MAX;
        min_idx = -1;
        for(j = 0; j < comm_sz; j++){
            if(merge_indexs[j] < merge_ends[j] && exchanged_data[merge_indexs[j]] <= min){
                min = exchanged_data[merge_indexs[j]];
                min_idx = j;
            }
        }
        (*merged_data)[i] = min;
        merge_indexs[min_idx]++;
    }

    free(merge_indexs);
    free(merge_ends);
    free(merge_counts);
    free(exchanged_data);
}

void gather_data(int sub_count, long* sub_data, long* data, int my_rank, int comm_sz){
    int i, *recv_counts = (int*)malloc(comm_sz * sizeof(int)), *recv_disps = (int*)malloc(comm_sz * sizeof(int));
    MPI_Gather(&sub_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(my_rank == 0){
        recv_disps[0] = 0;
        for(i = 1; i < comm_sz; i++){
            recv_disps[i] = recv_disps[i-1] + recv_counts[i-1];
        }
    }

    MPI_Gatherv(sub_data, sub_count, MPI_LONG, data, recv_counts, recv_disps, MPI_LONG, 0, MPI_COMM_WORLD);

    free(sub_data);
    free(recv_counts);
    free(recv_disps);
}

int main(int argc, char** argv){
    // 变量声明
    int my_rank, comm_sz, *merge_counts;
    size_t i = 0, init_data_count = 0, data_count = 0, local_data_count = 0, exchanged_data_count = 0;
    long *data = NULL, *local_data = NULL, *pivots = NULL, *exchanged_data = NULL, *merged_data = NULL;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // 数据读取
    if(my_rank == 0) {
        read_from_file("psrs_data", comm_sz, &init_data_count, &data_count, &data);
    }
    printf("Process %d read_from_file finished!\n", my_rank);

    // 数据分发
    deliver_data(data_count, data, my_rank, comm_sz, &local_data_count, &local_data);
    printf("Process %d deliver_data finished!\n", my_rank);

    // 本地排序+
    qsort(local_data, local_data_count, sizeof(long), cmp);
    printf("Process %d qsort finished!\n", my_rank);

    // 数据采样
    select_pivot(local_data_count, local_data, my_rank, comm_sz, &pivots);
    printf("Process %d select_pivot finished!\n", my_rank);

    // 数据交换
    exchange_data(local_data_count, local_data, pivots, comm_sz, &exchanged_data_count, &exchanged_data, &merge_counts);
    printf("Process %d exchange_data finished!\n", my_rank);

    // 数据归并
    merge_data(exchanged_data_count, exchanged_data, comm_sz, merge_counts,  &merged_data);
    printf("Process %d merge_data finished!\n", my_rank);

    // 数据汇集
    gather_data((int)exchanged_data_count, merged_data, data, my_rank, comm_sz);
    printf("Process %d gather_data finished!\n", my_rank);

    printf("Process %d all finished!\n", my_rank);

    if(my_rank == 0){
        for(i = 0; i < init_data_count; i++){
            if(i > 0 && data[i] < data[i-1]){
                printf("Error!\n");
                break;
            }
//            printf("%ld\n", data[i]);
        }
        free(data);
    }

    MPI_Finalize();
    return 0;
}
