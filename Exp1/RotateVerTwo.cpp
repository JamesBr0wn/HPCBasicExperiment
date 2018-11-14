#include <stdio.h>
#include <stdlib.h>

#define RIDX(row, col, dim) (row) * (dim) + col
#define COPY(dst, src) *(dst) = *(src)

typedef struct pixel{
  int red;
  int green;
  int blue;
}pixel;

void generate_graph(int dim, pixel** src, pixel** dst){
  int i, j;
  *src = (pixel*)malloc(sizeof(pixel) * dim * dim);
  *dst = (pixel*)malloc(sizeof(pixel) * dim * dim);
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      (*src)[RIDX(i, j, dim)].red = (int)rand();
      (*src)[RIDX(i, j, dim)].green = (int)rand();
      (*src)[RIDX(i, j, dim)].blue = (int)rand();
    }
  }
}

void free_space(pixel* src, pixel* dst){
  free(src);
  free(dst);
}

void show(int dim, pixel* graph){
  int i, j;
  printf("----------------------------------------\n");
  printf("Red:\n");
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      printf("%d ", graph[RIDX(i, j, dim)].red);
    }
    putchar('\n');
  }
  printf("----------------------------------------\n");
  printf("----------------------------------------\n");
  printf("Green:\n");
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      printf("%d ", graph[RIDX(i, j, dim)].green);
    }
    putchar('\n');
  }
  printf("----------------------------------------\n");
  printf("----------------------------------------\n");
  printf("Blue:\n");
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      printf("%d ", graph[RIDX(i, j, dim)].blue);
    }
    putchar('\n');
  }
  printf("----------------------------------------\n");
}

void rotate_ver_two(int dim, pixel* src, pixel* dst){
  int i, j, ii, jj;
  for(ii = 0; ii < dim; ii += 32){
    for(jj = 0; jj < dim; jj += 32){
      for(i = ii; i < ii + 32; i += 4){
        for(j = jj; j < jj + 32; j += 4){
          dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
          dst[RIDX(dim-1-j, i+1, dim)] = src[RIDX(i+1, j, dim)];
          dst[RIDX(dim-1-j, i+2, dim)] = src[RIDX(i+2, j, dim)];
          dst[RIDX(dim-1-j, i+3, dim)] = src[RIDX(i+3, j, dim)];
          dst[RIDX(dim-1-j-1, i, dim)] = src[RIDX(i, j+1, dim)];
          dst[RIDX(dim-1-j-1, i+1, dim)] = src[RIDX(i+1, j+1, dim)];
          dst[RIDX(dim-1-j-1, i+2, dim)] = src[RIDX(i+2, j+1, dim)];
          dst[RIDX(dim-1-j-1, i+3, dim)] = src[RIDX(i+3, j+1, dim)];
          dst[RIDX(dim-1-j-2, i, dim)] = src[RIDX(i, j+2, dim)];
          dst[RIDX(dim-1-j-2, i+1, dim)] = src[RIDX(i+1, j+2, dim)];
          dst[RIDX(dim-1-j-2, i+2, dim)] = src[RIDX(i+2, j+2, dim)];
          dst[RIDX(dim-1-j-2, i+3, dim)] = src[RIDX(i+3, j+2, dim)];
          dst[RIDX(dim-1-j-3, i, dim)] = src[RIDX(i, j+3, dim)];
          dst[RIDX(dim-1-j-3, i+1, dim)] = src[RIDX(i+1, j+3, dim)];
          dst[RIDX(dim-1-j-3, i+2, dim)] = src[RIDX(i+2, j+3, dim)];
          dst[RIDX(dim-1-j-3, i+3, dim)] = src[RIDX(i+3, j+3, dim)];
        }
      }
    }
  }
}

int main(){
  int dim = 4096;
  pixel *src, *dst;
  generate_graph(dim, &src, &dst);
  // show(dim, src);
  rotate_ver_two(dim, src, dst);
  // show(dim, dst);
  free_space(src, dst);
  return 0;
}
