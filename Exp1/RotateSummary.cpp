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

void naive_rotate(int dim, pixel* src, pixel* dst){
  int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
    }
  }
}

void rotate_ver_one(int dim, pixel* src, pixel* dst){
  int i, j, ii, jj;
  for(ii = 0; ii < dim; ii += 4){
    for(jj = 0; jj < dim; jj += 4){
      for(i = ii; i < ii + 4; i++){
        for(j = jj; j < jj + 4; j++){
          dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
        }
      }
    }
  }
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

void rotate_ver_three(int dim, pixel* src, pixel* dst){
  int i, j;
  for(i = 0; i < dim; i += 32){
    for(j = dim-1; j >= 0; j--){
      pixel* dptr = dst + RIDX(dim-1-j, i,  dim);
      pixel* sptr = src + RIDX(i, j, dim);
      COPY(dptr, sptr);
      sptr += dim;
      COPY(dptr+1, sptr);
      sptr += dim;
      COPY(dptr+2, sptr);
      sptr += dim;
      COPY(dptr+3, sptr);
      sptr += dim;
      COPY(dptr+4, sptr);
      sptr += dim;
      COPY(dptr+5, sptr);
      sptr += dim;
      COPY(dptr+6, sptr);
      sptr += dim;
      COPY(dptr+7, sptr);
      sptr += dim;
      COPY(dptr+8, sptr);
      sptr += dim;
      COPY(dptr+9, sptr);
      sptr += dim;
      COPY(dptr+10, sptr);
      sptr += dim;
      COPY(dptr+11, sptr);
      sptr += dim;
      COPY(dptr+12, sptr);
      sptr += dim;
      COPY(dptr+13, sptr);
      sptr += dim;
      COPY(dptr+14, sptr);
      sptr += dim;
      COPY(dptr+15, sptr);
      sptr += dim;
      COPY(dptr+16, sptr);
      sptr += dim;
      COPY(dptr+17, sptr);
      sptr += dim;
      COPY(dptr+18, sptr);
      sptr += dim;
      COPY(dptr+19, sptr);
      sptr += dim;
      COPY(dptr+20, sptr);
      sptr += dim;
      COPY(dptr+21, sptr);
      sptr += dim;
      COPY(dptr+22, sptr);
      sptr += dim;
      COPY(dptr+23, sptr);
      sptr += dim;
      COPY(dptr+24, sptr);
      sptr += dim;
      COPY(dptr+25, sptr);
      sptr += dim;
      COPY(dptr+26, sptr);
      sptr += dim;
      COPY(dptr+27, sptr);
      sptr += dim;
      COPY(dptr+28, sptr);
      sptr += dim;
      COPY(dptr+29, sptr);
      sptr += dim;
      COPY(dptr+30, sptr);
      sptr += dim;
      COPY(dptr+31, sptr);
      sptr += dim;
    }
  }
}

int main(){
  int dim = 4096;
  pixel *src, *dst;
  generate_graph(dim, &src, &dst);
  // show(dim, src);
  naive_rotate(dim, src, dst);
  free_space(src, dst);
  generate_graph(dim, &src, &dst);
  rotate_ver_one(dim, src, dst);
  free_space(src, dst);
  generate_graph(dim, &src, &dst);
  rotate_ver_two(dim, src, dst);
  free_space(src, dst);
  generate_graph(dim, &src, &dst);
  rotate_ver_three(dim, src, dst);
  // show(dim, dst);
  free_space(src, dst);
  return 0;
}
