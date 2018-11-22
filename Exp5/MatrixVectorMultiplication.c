#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int my_rank, comm_sz;

void get_size(int* m_ptr, int* local_m_ptr, int* n_ptr, int* local_n_ptr){
    if(my_rank == 0){
        printf("Please enter the size of matrix and vector:\t");
        scanf("%d", n_ptr);
    }
    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    *m_ptr = *n_ptr;
    *local_m_ptr = *m_ptr / comm_sz;
    *local_n_ptr = *n_ptr / comm_sz;
}

void allocate_space(int m, int local_m, int local_n, double** local_mat_p, double** local_vec_p,  double** local_res_p){
    *local_mat_p = (double*)malloc(m * local_n * sizeof(double));
    *local_vec_p = (double*)malloc(local_n * sizeof(double));
    *local_res_p = (double*)malloc(local_m * sizeof(double));
}

void create_type(int m, int n, int local_n, MPI_Datatype* local_mat_t_ptr){
    MPI_Datatype mpi_vec_t;
    MPI_Type_vector(m, local_n, n, MPI_DOUBLE, &mpi_vec_t);
    MPI_Type_create_resized(mpi_vec_t, 0, local_n * sizeof(double), local_mat_t_ptr);
    MPI_Type_commit(local_mat_t_ptr);
}

// Matrix shape: m * n
void read_matrix(int m, int n, int local_n, double *local_mat, MPI_Datatype local_mat_t){
    double* matrix = NULL;
    int i;

    if(my_rank == 0){
        matrix = (double*)malloc(m * n * sizeof(double));
        printf("Please enter the matrix:\n");
        for(i = 0; i < m * n; i++){
            scanf("%lf", matrix + i);
        }

        MPI_Scatter(matrix, 1, local_mat_t, local_mat, m * local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(matrix);
    }else{
        MPI_Scatter(matrix, 1, local_mat_t, local_mat, m * local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

// Vector shape: n * 1
void read_vector(int n, int local_n, double* local_vec){
    double* vector = NULL;
    int i;

    if(my_rank == 0){
        vector = (double*)malloc(n * sizeof(double));
        printf("Please enter the vector:\n");
        for(i  = 0; i < n; i++){
            scanf("%lf", vector + i);
        }

        MPI_Scatter(vector, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(vector);
    }else{
        MPI_Scatter(vector, local_n, MPI_DOUBLE, local_vec, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void data_compute(int m, int local_m, int local_n, double *local_mat, double *local_vec, double *local_res){
    int i, j, *recv_counts;
    double* median_res;

    median_res = (double*)malloc(m * sizeof(double));
    for(i = 0; i < m; i++){
        median_res[i] = 0;
        for(j = 0; j < local_n; j++){
            median_res[i] += local_mat[i * local_n + j] * local_vec[j];
        }
    }

    recv_counts = (int*)malloc(comm_sz * sizeof(int));
    for(i = 0; i < comm_sz; i++){
        recv_counts[i] = local_m;
    }

    MPI_Reduce_scatter(median_res, local_res, recv_counts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    free(median_res);
    free(recv_counts);
}

void print_matrix(char* mat_name, int m, int n, int local_n, double* local_mat, MPI_Datatype local_mat_t){
    double* matrix = NULL;
    int i, j;

    if(my_rank == 0){
        matrix = (double*)malloc(m * n * sizeof(double));

        MPI_Gather(local_mat, m * local_n, MPI_DOUBLE, matrix, 1, local_mat_t, 0, MPI_COMM_WORLD);

        printf("%s\n", mat_name);
        for(i = 0; i < m; i++){
            for(j = 0; j < n; j++){
                printf("%.3f ", matrix[i*n+j]);
            }
            putchar('\n');
        }

        free(matrix);
    }else{
        MPI_Gather(local_mat, m * local_n, MPI_DOUBLE, matrix, 1, local_mat_t, 0, MPI_COMM_WORLD);
    }
}

void print_vector(char*vec_name, int n, int local_n, double* local_vec){
    double* vector = NULL;
    int i;

    if(my_rank == 0){
        vector = (double*)malloc(n * sizeof(double));

        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        printf("%s\n", vec_name);
        for(i = 0; i < n; i++){
            printf("%.3lf ", vector[i]);
        }
        putchar('\n');

        free(vector);
    }else{
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv){
    int m, n, local_m, local_n;
    double *local_mat = NULL, *local_vec = NULL, *local_res = NULL;
    MPI_Datatype local_mat_t;

    // MPI初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    get_size(&m, &local_m, &n, &local_n);

    allocate_space(m, local_m, local_n, &local_mat, &local_vec, &local_res);

    create_type(m, n, local_n, &local_mat_t);

    read_matrix(m, n, local_n, local_mat, local_mat_t);

    print_matrix("A", m, n, local_n, local_mat, local_mat_t);

    read_vector(n, local_n, local_vec);

    print_vector("x", n, local_n, local_vec);

    data_compute(m, local_m, local_n, local_mat, local_vec, local_res);

    print_vector("y=Ax", m, local_m, local_res);

    free(local_mat);
    free(local_vec);
    free(local_res);

    MPI_Type_free(&local_mat_t);
    MPI_Finalize();

    return 0;
}
