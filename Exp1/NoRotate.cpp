#include <stdio.h>
#include <stdlib.h>

#define RIDX(row, col, dim) (row) * (dim) + col

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

int main(){
  int dim = 4096;
  pixel *src, *dst;
  generate_graph(dim, &src, &dst);
  // show(dim, src);
  // show(dim, dst);
  free_space(src, dst);
  return 0;
}
