//
// Created by chen on 3/14/2022.
//

#ifndef SOFTWARE_PROJECT_FINAL_KMEANS_H
#define SOFTWARE_PROJECT_FINAL_KMEANS_H

# include <Python.h>

double **initialize_2d_double_array(double *arr[], int d1, int d2);

int calc_argmin(double *mu[], double *vector, int k, int dimension);

double compute_distance(double vec1[], double vec2[], int dimension);

void reset_clusters(double **vectors_list, double *mu[], double *new_sum[], int k, int N, int dimension);

double calculating_epsilon(double *mu[], double *new_mu[], int k, int dimension);

int check_allocation(const double *p);

void free_memory(double **array, int len);

double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns);

PyObject *transform_2dArray_to_PyObject(double **mat, int rows, int columns);



#endif //SOFTWARE_PROJECT_FINAL_KMEANS_H
