#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

float* generateData(int data_count, float min_meas, float max_meas){
  // 产生模拟数据
  srand((unsigned int)(time(NULL)));
  float* data = (float*)malloc(data_count * sizeof(float));
  int i;
  for(i = 0; i < data_count; i++){
    data[i] = min_meas + (max_meas - min_meas) * rand() / (double)RAND_MAX;
  }
  return data;
}

float* calculateBinMaxes(float min_meas, float max_meas, int bin_count){
  // 计算直方图各桶数据范围
  int i;
  float bin_width = (max_meas - min_meas) / bin_count, current_max = min_meas + bin_width; 
  float* bin_maxes = (float*)malloc(bin_count * sizeof(float));
  for(i = 0; i < bin_count; i++){
    bin_maxes[i] = current_max;
    current_max += bin_width;
  }
  bin_maxes[bin_count-1] = max_meas + 1;
  return bin_maxes;
}

int* calculateHistogram(int my_id, int comm_sz, int data_count, float* data, int bin_count, float* bin_maxes){
  int i, j;
  int* histogram = (int*)malloc(bin_count * sizeof(int));
  
  // 计算自己的部分直方图
  memset(histogram, 0, bin_count * sizeof(int));
  for(i = 0; i < data_count; i++){
    for(j = 0; j < bin_count; j++){
      if(data[i] < bin_maxes[j]){
        histogram[j]++;
        break;
      }
    }
  }

  int comm_count = 1;   // 与通信轮数有关的计数值，用于判断一个进程是发送进程还是接受进程
  int* temp_hist = (int*)malloc(bin_count * sizeof(int));   // 存放接收到的部分直方图
  MPI_Status status;
  while(comm_count != comm_sz){
    comm_count *= 2;    // 计数值乘2，进入下一轮通信
    if(my_id % comm_count != 0){    // 对于发送进程，将自己对应的部分直方图发送出去，释放申请的内存空间并返回NULL
      MPI_Send(histogram, bin_count, MPI_INT, my_id / comm_count * comm_count, comm_count, MPI_COMM_WORLD);
      printf("Process %d send local histogram to %d\n", my_id, my_id / comm_count * comm_count);
      free(histogram);
      printf("Work of process %d has finished!\n", my_id);
      free(temp_hist);
      return NULL;
    }else{    // 对于接收进程，接收来自其他进程的部分直方图，与自己的部分直方图累加
      MPI_Recv(temp_hist, bin_count, MPI_INT, MPI_ANY_SOURCE, comm_count, MPI_COMM_WORLD, &status);
      printf("Process %d receieve local histogram from %d\n", my_id, status.MPI_SOURCE);
      for(j = 0; j < bin_count; j++){
        histogram[j] += temp_hist[j];
      }
    }
  }
  // 最后只有一个进程（进程0）得到完整的直方图，由它返回该完整直方图
  free(temp_hist);
  printf("Work of process %d has finished!\n", my_id);
  return histogram;
}

int main(int argc, char** argv){
  // 变量初始化
  int data_count = 100000000, bin_count = 10, my_id, comm_sz, i, local_data_count;
  int *temp_hist, *histogram;
  float min_meas = 0, max_meas = 10;
  float *local_data, *bin_maxes;

  // MPI初始化
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  bin_maxes = calculateBinMaxes(min_meas, max_meas, bin_count);
  if(my_id == 0){   // 对于进程0
    float* data = generateData(data_count, min_meas, max_meas);
    int start, count;
    int ele_per_pro = data_count / comm_sz;
    for(i = 0; i < comm_sz; i++){
      // 部分数据开始索引
      start = i * ele_per_pro;
      // 部分数据大小
      count = start + ele_per_pro < data_count ? ele_per_pro : data_count - start;
      if(i == 0){   // 对于自己的数据，直接赋值
        local_data_count = count;
        local_data = data;
      }else{    // 对于其他进程的数据，通过MPI_Send发送出去
        MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        printf("Process %d send local data count to %d\n", my_id, i);
        MPI_Send(data + start, count, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        printf("Process %d send local data to %d\n", my_id, i);
      }
    }
  }else{    // 对于其他进程
    MPI_Status status;
    // 接收部分数据大小
    MPI_Recv(&local_data_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    printf("Process %d receieve local data count from %d\n", my_id, status.MPI_SOURCE);
    // 根据部分数据大小申请存储位置
    local_data = (float*)malloc(local_data_count * sizeof(float));
    // 接收部分数据
    MPI_Recv(local_data, local_data_count, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
    printf("Process %d receieve local data from %d\n", my_id, status.MPI_SOURCE);
  }

  // 统计直方图
  temp_hist = calculateHistogram(my_id, comm_sz, local_data_count, local_data, bin_count, bin_maxes);

  // 同步
  MPI_Barrier(MPI_COMM_WORLD);

  // 结果打印和资源释放
  if(my_id == 0){
    float bin_start, bin_end;
    histogram = temp_hist;
    for(i = 0; i < bin_count; i++){
      if(i == 0){
        bin_start = min_meas;
      }else{
        bin_start = bin_maxes[i-1];
      }
      if(i == bin_count - 1){
        bin_end = max_meas;
      }else{
        bin_end = bin_maxes[i];
      }
      printf("%.2f~%.2f:\t%d\n", bin_start, bin_end, histogram[i]);
    }
    free(histogram);
    free(bin_maxes);
    free(local_data);
  }else{
    free(temp_hist);
    free(local_data);
  }
  MPI_Finalize();
  return 0;
}
