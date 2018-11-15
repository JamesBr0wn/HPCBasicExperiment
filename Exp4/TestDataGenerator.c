#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
        printf("%ld\n", temp);
        sprintf(str, "%ld\n", temp);
        fwrite(&temp, sizeof(long), 1, output_bin);
        fwrite(&str, strlen(str), 1, output_str);
    }
    fclose(output_bin);
    fclose(output_str);
    printf("Finished!\n");
    return 0;
}