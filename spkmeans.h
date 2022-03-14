//
// Created by chen on 3/12/2022.
//

#ifndef SOFTWARE_PROJECT_FINAL_SPKMEANS_H
#define SOFTWARE_PROJECT_FINAL_SPKMEANS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
double **sp_kmeans(double **data_points, int n, int dimension, int k);
double **weight_matrix(double **data_points, int n, int dimension);
double **diagonal_d_matrix(double **data_points, int n, int dimension);
double **laplacian_Lnorm(double **data_points, int n, int dimension);
double **jacobi_algo(double **data_points, int n, int dimension);
#endif //SOFTWARE_PROJECT_FINAL_SPKMEANS_H
