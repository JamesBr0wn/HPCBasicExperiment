#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>

int main(){
    FILE* file;
    int n, i, j;
    double temp;

    printf("Please enter the size of matrix/vector to be generated:\t");
    scanf("%d", &n);

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

    return 0;
}