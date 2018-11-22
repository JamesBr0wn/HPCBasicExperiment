#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

int cmp ( const void *a , const void *b ) {
    return *(long*)a > *(long*)b ? 1 : -1;
}

int main(){
    FILE* output_bin = fopen("psrs_data", "wb");
    FILE* output_str = fopen("psrs_str", "w");
    long temp, i, n;
    char str[20];
    printf("Please enter the amount of number to be generated:\t");
    scanf("%ld", &n);
    fwrite(&n, sizeof(long), 1, output_bin);
    sprintf(str, "%ld\n", n);
    fwrite(&str, strlen(str), 1, output_str);
    for(i = 0; i < n; i++){
        temp = random();
        sprintf(str, "%ld\n", temp);
        fwrite(&temp, sizeof(long), 1, output_bin);
        fwrite(&str, strlen(str), 1, output_str);
    }
    fclose(output_bin);
    fclose(output_str);
    printf("Finished!\n");

//    size_t init_data_count, data_count;
//    long* data;
//    int comm_sz = 4;
//    FILE* input_file = fopen("psrs_data", "rb");
//    fread(&init_data_count, sizeof(long), 1, input_file);
//
//    data_count = init_data_count % comm_sz == 0 ? init_data_count : init_data_count  + comm_sz - init_data_count % comm_sz;
//
//
//    data = (long*)malloc(init_data_count * sizeof(long));
//    fread(data, sizeof(long), init_data_count, input_file);
//
//    for(i = data_count - 1; i >= init_data_count; i--){
//        data[i] = LONG_MAX;
//    }
//    fclose(input_file);
//
//    qsort(data, data_count, sizeof(long), cmp);
//    for(i = 0; i < data_count; i++){
//        printf("%ld\n", data[i]);
//    }
    return 0;
}
