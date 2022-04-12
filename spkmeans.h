//
// Created by chen on 3/12/2022.
//

#ifndef SOFTWARE_PROJECT_FINAL_SPKMEANS_H
#define SOFTWARE_PROJECT_FINAL_SPKMEANS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void free_mat(double **matrix, int n);
void print_matrix(double **matrix, int rows, int cols);
double **weight_matrix(double **data_points, int n, int dimension);
double **diagonal_d_matrix(double **data_points, int n, int dimension);
double **laplacian_Lnorm(double **data_points, int n, int dimension);
int eigengap_heuristic(double *eigen_vals, int n);
double **jacobi_algo(double **data_points, int n);
double* transpose_and_sort(double **jacobi_mat, int n);
double **create_T_matrix(double **matrix, int k, int n);

#endif //SOFTWARE_PROJECT_FINAL_SPKMEANS_H
