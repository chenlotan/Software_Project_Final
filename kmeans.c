#include <math.h>
#include "kmeans.h"


double compute_distance(double vec1[], double vec2[], int dimension) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return sum;
}


double **initialize_2d_double_array(double **arr, int d1, int d2) {
    int i;
    arr = (double **) malloc(d1 * sizeof(double *));
    for (i = 0; i < d1; ++i) {
        arr[i] = (double *) calloc(d2, sizeof(double));
    }
    return arr;
}

int calc_argmin(double *mu[], double *vector, int k, int dimension) {
    double min_val = HUGE_VAL, sum_p;
    int min_mu = 0, i;
    /*printf("vector \n");
    for (i=0; i<dimension; ++i){
        printf("%f ,", vector[i]);
    }*/

    for (i = 0; i < k; ++i) {
        sum_p = compute_distance(mu[i], vector, dimension);
        if (sum_p < min_val - 0.000001) {
            min_val = sum_p;
            min_mu = i;
            /*printf("i = %d, minVal = %f , sump = %f \n", i, min_val, sum_p);*/

        }
    }
    return min_mu;
}

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[], int k, int N, int dimension) {
    int *count;
    double *vec;
    int t, j, r, s, q, min_mu;
    count = calloc(k, sizeof(int));

    for (t = 0; t < N; t++) {

        min_mu = calc_argmin(mu, vectors_list[t], k, dimension);
        count[min_mu]++;
        vec = vectors_list[t];
        for (j = 0; j < dimension; ++j) {
            new_sum[min_mu][j] += *vec;
            vec++;
        }
        /*printf("new_sum\n");
        for (int i=0; i<k; ++i){
            for (j=0; j<dimension; ++j){
                printf("%f, ", new_sum[i][j]);
            }
            printf("\n");
        }*/
    }
    for (r = 0; r < k; ++r) {
        if (count[r] == 0) {
            for (s = 0; s < dimension; ++s) {
                new_sum[r][s] = mu[r][s];
            }
        } else {
            for (q = 0; q < dimension; ++q) {
                new_sum[r][q] = new_sum[r][q] / (double) count[r];
            }
        }
    }
    free(count);
}

double calculating_epsilon(double *mu[], double *new_mu[], int k, int dimension) {
    double eps = 0, dist;
    int i;
    for (i = 0; i < k; ++i) {
        dist = compute_distance(mu[i], new_mu[i], dimension);
        if (eps < dist) {
            eps = dist;
        }
    }
    return eps;
}


int check_allocation(const double *p) {
    if (p == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    return 0;
}


void free_memory(double **array, int len) {
    int q;
    for (q = 0; q < len; q++) {
        free(array[q]);
    }
    free(array);
}




