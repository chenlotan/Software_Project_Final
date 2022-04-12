
#include "spkmeans.h"
#include "kmeans.h"

static PyObject *spk_ext(PyObject *self, PyObject *args);
static PyObject *fit(PyObject *self, PyObject *args);

static PyMethodDef k_means_func[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, NULL},
        {"spk_ext", (PyCFunction) spk_ext, METH_VARARGS, NULL},
        {NULL, NULL, 0,                          NULL}
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        k_means_func
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m)
        return NULL;
    return m;
}

static PyObject *fit(PyObject *self, PyObject *args) {
    int max_iter, i, j, q, k, dimension, N;
    double eps, curr_eps, **centroids, **data_points, **final_centroids = NULL;
    PyObject *centroids_copy, *data_points_copy;
    if (!PyArg_ParseTuple(args, "iiiidOO", &k, &dimension, &N, &max_iter, &eps, &centroids_copy, &data_points_copy)) {
        return NULL;
    }
    if (!(PyList_Check(centroids_copy)) || !(PyList_Check(data_points_copy))) {
        printf("Error in Types of Arguments\n");
    }
    final_centroids = initialize_2d_double_array(final_centroids, k, dimension);
    centroids = transform_PyObject_to_2dArray(centroids_copy, k, dimension);
    data_points = transform_PyObject_to_2dArray(data_points_copy, N, dimension);

    for (i = 0; i < max_iter; ++i) {
        reset_clusters(data_points, centroids, final_centroids, k, N, dimension);
        curr_eps = calculating_epsilon(centroids, final_centroids, k , dimension);
        for (j = 0; j < k; ++j) {
            for (q = 0; q < dimension; ++q) {
                centroids[j][q] = final_centroids[j][q];
                final_centroids[j][q] = 0;
            }
        }
        if (curr_eps < 0.000001) {
            break;
        }
    }
    centroids_copy = transform_2dArray_to_PyObject(centroids, k, dimension);
    free_memory(final_centroids, k);
    free_memory(centroids, k);
    free_memory(data_points, N);
    return centroids_copy;
}

static PyObject *spk_ext(PyObject *self, PyObject *args) {
    int k, goal, n, dimension;
    double **data_points, **result, **eigen_mat, **lnorm_mat, *eigen_vals;
    PyObject *data_points_copy, *result_py;
    if (!PyArg_ParseTuple(args, "iiiiO", &k, &goal, &n, &dimension, &data_points_copy)) {
        return NULL;
    }
    data_points = transform_PyObject_to_2dArray(data_points_copy, n, dimension);
    switch (goal) {
        case 1:
            lnorm_mat = laplacian_Lnorm(data_points, n, dimension);
            eigen_mat = jacobi_algo(lnorm_mat, n);
            eigen_vals = transpose_and_sort(eigen_mat, n);
            if (k==0){
                k = eigengap_heuristic(eigen_vals, n);
//                if (k==1){
//                    printf("An Error Has Occurred");
//                    exit(1);
//                }
            }
            result = create_T_matrix(eigen_mat, k, n);
            result_py = transform_2dArray_to_PyObject(result, n ,k);
            free_memory(lnorm_mat, n);
            free_memory(eigen_mat, n+1);
            free_memory(result, k);
            free_memory(data_points,n);
            break;
        case 2:
            result = weight_matrix(data_points, n, dimension);
            result_py = transform_2dArray_to_PyObject(result, n ,n);
            free_memory(result, n);
            free_memory(data_points,n);
            break;
        case 3:
            result = diagonal_d_matrix(data_points, n, dimension);
            result_py = transform_2dArray_to_PyObject(result, n ,n);
            free_memory(result, n);
            free_memory(data_points,n);
            break;
        case 4:
            result = laplacian_Lnorm(data_points, n, dimension);
            result_py = transform_2dArray_to_PyObject(result, n ,n);
            free_memory(result, n);
            free_memory(data_points,n);
            break;
        case 5:
            result = jacobi_algo(data_points, n);
            result_py = transform_2dArray_to_PyObject(result, n ,n);
            free_memory(result, n);
            free_memory(data_points,n);
            break;
        default:
            printf("Invalid Input!");
            exit(1);
    }

    return result_py;
}

double **transform_PyObject_to_2dArray(PyObject *mat, int rows, int columns) {
    double **new_mat = NULL;
    PyObject *row, *column;
    int i, j;
    new_mat = initialize_2d_double_array(new_mat, rows, columns);
    for (i = 0; i < rows; ++i) {
        row = PyList_GetItem(mat, i);
        for (j = 0; j < columns; ++j) {
            column = PyList_GetItem(row, j);
            new_mat[i][j] = PyFloat_AsDouble(column);
        }
    }
    return new_mat;
}

PyObject *transform_2dArray_to_PyObject(double **mat, int rows, int columns) {
    PyObject *new_mat, *row;
    int i, j;
    new_mat = PyList_New(rows);
    for (i = 0; i < rows; ++i) {
        row = PyList_New(columns);
        for (j = 0; j < columns; ++j) {
            PyList_SetItem(row, j, Py_BuildValue("f", mat[i][j]));
        }
        PyList_SetItem(new_mat, i, row);
    }
    return new_mat;
}

static PyMethodDef Kmeans_Methods[] = {
        {"fit",
                (PyCFunction) fit,
                METH_VARARGS,
                        NULL},
        {"spk_ext",
            (PyCFunction) spk_ext,
                     METH_VARARGS,
                         NULL},
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        Kmeans_Methods
};

PyMODINIT_FUNC
PyInit_kmeans(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
