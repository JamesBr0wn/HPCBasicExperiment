#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int my_rank, my_row_rank, my_col_rank, my_pivot_rank, comm_sz, row_comm_sz, col_comm_sz, pivot_comm_sz;
MPI_Comm total_comm, row_comm, col_comm, pivot_comm;

// 创建通信域
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

    // 创建笛卡尔拓扑通
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &total_comm);

    // 创建行通信域(保持行号)
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(total_comm, remain_dims, &row_comm);
    MPI_Comm_size(row_comm, &row_comm_sz);
    MPI_Comm_rank(row_comm, &my_row_rank);

    // 创建列通信域（保持列号）
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(total_comm, remain_dims, &col_comm);
    MPI_Comm_size(col_comm, &col_comm_sz);
    MPI_Comm_rank(col_comm, &my_col_rank);

    // 创建主元通信域
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
    // 进程0读入矩阵和向量大小
    if(my_rank == 0){
        printf("Please enter the size of matrix and vector:\t");
        scanf("%d", n_ptr);
    }

    // 将大小广播到所有进程
    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 计算线程局部矩阵、向量大小
    *m_ptr = *n_ptr;
    *local_m_ptr = *m_ptr / row_comm_sz;
    *local_n_ptr = *n_ptr / col_comm_sz;
}

// 产生模拟数据
void generate_data(int n){
    if(my_rank == 0){
        FILE* file;
        int i, j;
        double temp;
        file = fopen("matrix.txt", "w");
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                temp = rand() / (double)RAND_MAX * 1024;
                fprintf(file, "%lf ", temp);
            }
            fputc('\n', file);
        }
        fclose(file);

        file = fopen("vector.txt", "w");
        for(i = 0; i < n; i++){
            temp = rand() / (double)RAND_MAX * 1024;
            fprintf(file, "%lf ", temp);
        }
        fputc('\n', file);
        fclose(file);
    }
}

// 为各进程局部运算数和结果分配内存空间
void allocate_space(int local_m, int local_n, double** local_mat_ptr, double** local_vec_ptr, double** local_res_ptr){
    *local_mat_ptr = (double*)malloc(local_m * local_n * sizeof(double));
    *local_vec_ptr = (double*)malloc(local_n * sizeof(double));
    *local_res_ptr = (double*)malloc(local_m * sizeof(double));
}

// 为按列划分矩阵创建自定义类型
void create_type(int local_m, int n, int local_n, MPI_Datatype* local_mat_t_ptr){
    MPI_Datatype mpi_vec_t;
    MPI_Type_vector(local_m, local_n, n, MPI_DOUBLE, &mpi_vec_t);
    MPI_Type_create_resized(mpi_vec_t, 0, local_n * sizeof(double), local_mat_t_ptr);
    MPI_Type_commit(local_mat_t_ptr);
}

// 读入矩阵
void read_matrix(int m, int local_m, int n, int local_n, double* local_mat, MPI_Datatype local_mat_t){
    double *matrix = NULL, *median_mat = NULL;
    int i;

    if(my_pivot_rank == 0) {    // 进程0
        // 读入矩阵
        matrix = (double *) malloc(m * n * sizeof(double));
//        printf("Please enter the matrix:\n");
//        for (i = 0; i < m * n; i++) {
//            scanf("%lf", matrix + i);
//        }
        FILE* file = fopen("matrix.txt", "r");
        for(i = 0; i < m * n; i++){
            fscanf(file, "%lf", matrix + i);
        }

        // 按行将矩阵划分为中间矩阵（按行划分），散发给每一行的主元进程
        median_mat = (double*)malloc(local_m * n * sizeof(double));
        MPI_Scatter(matrix, local_m * n, MPI_DOUBLE, median_mat, local_m * n, MPI_DOUBLE, 0, pivot_comm);
        free(matrix);

        // 按列将中间矩阵划分为部分矩阵（按块划分），散发给同一行的每个进程
        MPI_Scatter(median_mat, 1, local_mat_t, local_mat, local_m * local_n, MPI_DOUBLE, my_row_rank, col_comm);
        free(median_mat);
    }else if(my_pivot_rank != -1){  // 其他主元进程
        // 从进程0接收散发的中间矩阵
        median_mat = (double*)malloc(local_m * n * sizeof(double));
        MPI_Scatter(matrix, local_m * n, MPI_DOUBLE, median_mat, local_m * n, MPI_DOUBLE, 0, pivot_comm);

        // 按列将中间矩阵划分为部分矩阵（按块划分），散发给同一行的每个进程
        MPI_Scatter(median_mat, 1, local_mat_t, local_mat, local_m * local_n, MPI_DOUBLE, my_row_rank, col_comm);
        free(median_mat);
    }else{  // 其他非主元进程
        // 接收自己的部分矩阵
        MPI_Scatter(median_mat, 1, local_mat_t, local_mat, local_m * local_n, MPI_DOUBLE, my_row_rank, col_comm);
    }
}

// 打印矩阵
void print_matrix(char* mat_name, int m, int local_m, int n, int local_n, double* local_mat, MPI_Datatype local_mat_t){
    double *matrix = NULL, *median_mat = NULL;
    int i, j;

    if(my_pivot_rank == 0){ // 进程0
        // 一次聚集得到中间矩阵
        median_mat = (double*)malloc(local_m * n * sizeof(double));
        MPI_Gather(local_mat, local_m * local_n, MPI_DOUBLE, median_mat, 1, local_mat_t, my_row_rank, col_comm);

        // 二次聚集得到完整矩阵
        matrix = (double*)malloc(m * n * sizeof(double));
        MPI_Gather(median_mat, local_m * n, MPI_DOUBLE, matrix, local_m * n, MPI_DOUBLE, 0, pivot_comm);

        // 打印矩阵
        printf("%s\n", mat_name);
        for(i = 0; i < m; i++){
            for(j = 0; j < n; j++){
                printf("%.3f ", matrix[i*n+j]);
            }
            putchar('\n');
        }

        free(median_mat);
        free(matrix);
    }else if(my_pivot_rank != -1){  // 其他主元进程
        // 一次聚集得到中间矩阵
        median_mat = (double*)malloc(local_m * n * sizeof(double));
        MPI_Gather(local_mat, local_m * local_n, MPI_DOUBLE, median_mat, 1, local_mat_t, my_row_rank, col_comm);

        // 二次聚集得到完整矩阵
        MPI_Gather(median_mat, local_m * n, MPI_DOUBLE, matrix, local_m * n, MPI_DOUBLE, 0, pivot_comm);

        free(median_mat);
    }else{  // 其他非主元进程
        // 一次聚集得到中间矩阵
        MPI_Gather(local_mat, local_m * local_n, MPI_DOUBLE, median_mat, 1, local_mat_t, my_row_rank, col_comm);
    }
}

// 读入向量
void read_vector(int n, int local_n, double* local_vec){
    double* vector = NULL;
    int i;

    if(my_pivot_rank == 0) {    // 进程0
        // 读入向量
        vector = (double*)malloc(n * sizeof(double));
//        printf("Please enter the vector:\n");
//        for(i  = 0; i < n; i++){
//            scanf("%lf", vector + i);
//        }
        FILE* file = fopen("vector.txt", "r");
        for(i = 0; i < n; i++){
            fscanf(file, "%lf", vector + i);
        }

        // 向所有主元进程广播向量
        MPI_Bcast(vector, n, MPI_DOUBLE, 0, pivot_comm);

        // 作为主元进程向本行所有进程散发部分向量
        MPI_Scatter(vector, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, my_row_rank, col_comm);

        free(vector);
    }else if(my_pivot_rank != -1){  // 其他主元进程
        // 接收进程0广播的向量
        vector = (double*)malloc(n * sizeof(double));
        MPI_Bcast(vector, n, MPI_DOUBLE, 0, pivot_comm);

        // 作为主元进程向本行所有进程散发部分向量
        MPI_Scatter(vector, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, my_row_rank, col_comm);

        free(vector);
    }else{  // 其他非主元进程
        // 接收主元进程散发的部分向量
        MPI_Scatter(vector, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, my_row_rank, col_comm);
    }
}

// 打印向量
void print_vector(char*vec_name, int n, int local_n, double* local_vec){
    double* vector = NULL;
    int i;

    if(my_pivot_rank == 0){ // 进程0
        // 聚集得到完整向量
        vector = (double*)malloc(n * sizeof(double));
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, my_row_rank, col_comm);

        // 打印向量
        printf("%s\n", vec_name);
        for(i = 0; i < n; i++){
            printf("%.3lf ", vector[i]);
        }
        putchar('\n');

        free(vector);
    }else if(my_pivot_rank != -1){  // 其他主元进程
        // 聚集得到完整向量
        vector = (double*)malloc(n * sizeof(double));
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, my_row_rank, col_comm);

        free(vector);
    }else{  // 其他非主元进程
        // 向主元进程聚集部分向量
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, my_row_rank, col_comm);
    }
}

// 矩阵向量乘
void data_compute(int local_m, int local_n, double *local_mat, double *local_vec, double *local_res){
    int i, j;
    double* median_res;

    // 每个进程计算中间结果向量
    median_res = (double*)malloc(local_m * sizeof(double));
    for(i = 0; i < local_m; i++){
        median_res[i] = 0;
        for(j = 0; j < local_n; j++){
            median_res[i] += local_mat[i * local_n + j] * local_vec[j];
        }
    }

    // 聚集到主元进程得到部分结果向量
    MPI_Reduce(median_res, local_res, local_m, MPI_DOUBLE, MPI_SUM, my_row_rank, col_comm);

    free(median_res);
}

// 打印结果
void print_result(char* res_name, int m, int local_m, double* local_res){
    double* result = NULL;
    int i;

    if(my_pivot_rank == 0){ // 进程0
        // 聚集得到完整结果向量
        result = (double*)malloc(m * sizeof(double));
        MPI_Gather(local_res, local_m, MPI_DOUBLE, result, local_m, MPI_DOUBLE, 0, pivot_comm);

        // 打印结果向量
        printf("%s\n", res_name);
        for(i = 0; i < m; i++){
            printf("%.3lf ", result[i]);
        }
        putchar('\n');

        free(result);
    }else if(my_pivot_rank != -1){
        // 向进程0聚集以得到完整结果向量
        MPI_Gather(local_res, local_m, MPI_DOUBLE, result, local_m, MPI_DOUBLE, 0, pivot_comm);
    }
}

int main(int argc, char** argv){
    int m, n, local_m, local_n;
    double *local_mat, *local_vec, *local_res;
    MPI_Datatype local_mat_t;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    create_comm();

    get_size(&m, &local_m, &n, &local_n);

    generate_data(n);

    allocate_space(local_m, local_n, &local_mat, &local_vec, &local_res);

    create_type(local_m, n, local_n, &local_mat_t);

    read_matrix(m, local_m, n, local_n, local_mat, local_mat_t);

    print_matrix("A", m, local_m, m, local_n, local_mat, local_mat_t);

    read_vector(n, local_n, local_vec);

    print_vector("x", n, local_n, local_vec);

    data_compute(local_m, local_n, local_mat, local_vec, local_res);

    print_result("y", m, local_m, local_res);

    free(local_mat);
    free(local_vec);
    free(local_res);

    MPI_Type_free(&local_mat_t);
    MPI_Finalize();

    return 0;
}