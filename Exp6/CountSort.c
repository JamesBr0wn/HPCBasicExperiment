#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void Count_sort(int a[], int n) {
    int i, j, count;
    int* temp = malloc(n*sizeof(int));

    for (i = 0; i < n; i++) {
        count = 0;
        for (j = 0; j < n; j++)
            if (a[j] < a[i])
                count++;
            else if (a[j] == a[i] && j < i)
                count++;
        temp[count] = a[i];
    }

    memcpy(a, temp, n*sizeof(int));
    free(temp);
}

void Count_sort_omp(int a[], int n) {
    int i, j, count;
    int *temp = malloc(n * sizeof(int));
    #pragma omp parallel private(i, j, count) shared(n, a, temp)
    {
        #pragma omp for
        for (i = 0; i < n; i++) {
            count = 0;
            for (j = 0; j < n; j++)
                if (a[j] < a[i])
                    count++;
                else if (a[j] == a[i] && j < i)
                    count++;
            temp[count] = a[i];
        }
        #pragma omp for
        for(i = 0; i < n; i++){
            a[i] = temp[i];
        }
    }
    free(temp);
}

int main(){
    int i;
    int* a = malloc(10000 * sizeof(int));

    omp_set_num_threads(4);

    for(i = 0; i < 10000; i++){
        a[i] = rand();
    }

    for(i = 0; i < 20; i++){
        printf("%d\n", a[i]);
    }
    putchar('\n');


    Count_sort_omp(a, 10000);

    for(i = 0; i < 20; i++){
        printf("%d\n", a[i]);
    }
    putchar('\n');

    return 0;
}