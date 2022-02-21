#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int dimension, n;
double **data_points;


int main() {
    return 0;
}
/**
 * 1. check allocation
 * 2.
 *
 * **/


/** free memory of nxn matrix **/
void free_mat(double **matrix) {
    int i;
    for (i = 0; i < n; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

/** returns pointer to a matrix that initialized to zeros **/
double **generate_matrix() {
    int i;
    double **matrix;
    matrix = (double **) calloc(n, sizeof(double *));
    for (i = 0; i < n; ++i) {
        matrix[i] = (double *) calloc(n, sizeof(double));
    }
    return matrix;
}

/** sets the given matrix to an identity matrix **/
void set_identity(double **mat) {
    int i, j;
    for (i = 0; i < n; ++i) {
        mat[i][i] = 1;
        for (j = i + 1; j < n; ++j) {
            mat[i][j] = mat[j][i] = 0;
        }
    }
}

/** given two initialized matrices at size n, copying the data from the first matrix to the second one **/
void copy_matrix(double **matrix, double **copy_matrix) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            copy_matrix[i][j] = matrix[i][j];
        }
    }
}

/** multiply matrix V in matrix P in place **/
void multiply(double **V, double **P) {
    int i, j, k;
    double **A;
    A = generate_matrix();
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k) {
                A[i][j] += (V[i][k] * P[k][j]);
            }
        }
    }
    copy_matrix(A, V);
    free_mat(A);
}

/** computes the Euclidean norm between two vectors **/
double compute_weight(double vec1[], double vec2[]) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return exp(-sqrt(sum) / 2);
}

/** calculate the weight matrix out of the data points **/
double **weight_matrix() {
    int i, j;
    double **w_mat = generate_matrix();
    for (i = 0; i < n; ++i) {
        w_mat[i][i] = 0;
        for (j = i + 1; j < n; ++j) {
            w_mat[i][j] = compute_weight(data_points[i], data_points[j]);
            w_mat[j][i] = w_mat[i][j];
        }
    }
    return w_mat;
}

/** computes the diagonal matrix out of the weight matrix **/
double **diagonal_d_matrix() {
    double **diag_mat, **w_mat, sum;
    int i, j;
    diag_mat = generate_matrix();
    w_mat = weight_matrix();
    for (i = 0; i < n; ++i) {
        sum = 0;
        for (j = 0; j < n; ++j) {
            sum += w_mat[i][j];
        }
        diag_mat[i][i] = sum;
    }
    return diag_mat;
}

/** calculate Lnorm matrix - L_norm = I - D^(-0.5)WD^(-0.5) **/
double **laplacian_Lnorm() {
    int i, j;
    double **diag_mat, **w_mat, **Lnorm_mat;

    w_mat = weight_matrix();
    diag_mat = diagonal_d_matrix();
    Lnorm_mat = generate_matrix();

    for (i = 0; i < n; ++i) {
        Lnorm_mat[i][i] = 1 - (w_mat[i][i] / diag_mat[i][i]);
        for (j = i + 1; j < n; ++j) {
            Lnorm_mat[i][j] = w_mat[i][j] / sqrt(diag_mat[i][i] * diag_mat[j][j]);
            Lnorm_mat[j][i] = Lnorm_mat[i][j];
        }
    }
    free_mat(diag_mat);
    free_mat(w_mat);
    return Lnorm_mat;
}

/** calculate the f(A)^2 of a given matrix **/
double get_of_f_A(double **A) {
    int i, j;
    double sum = 0;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            sum += pow(A[i][j], 2);
        }
    }
    return 2 * sum;
}

/** calculate A' out of A, c and s **/
void cal_newA(double **new_A, double **A, double c, double s, int max_i, int max_j) {
    int r;
    for (r = 0; r < n; ++r) {
        if (r != max_i && r != max_j) {
            new_A[r][max_i] = c * A[r][max_i] - s * A[r][max_j];
            new_A[r][max_j] = c * A[r][max_j] + s * A[r][max_i];
        }
    }
    new_A[max_i][max_i] = pow(c, 2) * A[max_i][max_i] + pow(s, 2) * A[max_j][max_j] - 2 * s * c * A[max_i][max_j];
    new_A[max_j][max_j] = pow(c, 2) * A[max_j][max_j] + pow(s, 2) * A[max_i][max_i] + 2 * s * c * A[max_i][max_j];
    new_A[max_i][max_j] = new_A[max_j][max_i] = 0;
}

/** calculate c and s from the Lnorm matrix and returns 0 if it was valid and -1 elsewhere **/
int set_data(double **Lnorm_mat, double *c, double *s, int *max_i, int *max_j) {
    int i, j;
    double max_Aij = 0, theta, t;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            if (fabs(Lnorm_mat[i][j]) > max_Aij) {
                max_Aij = Lnorm_mat[i][j];
                *max_i = i;
                *max_j = j;
            }
        }
    }
    if (max_Aij == 0) {
        return -1;
    }
    theta = (Lnorm_mat[*max_j][*max_j] - Lnorm_mat[*max_i][*max_i]) / (2 * Lnorm_mat[*max_i][*max_j]);
    t = ((theta >= 0) - (theta < 0)) / (fabs(theta) + sqrt(pow(theta, 2)) + 1);
    *c = 1 / sqrt(pow(t, 2) + 1);
    *s = t * (*c);
    return 0;
}

/** Jacobi Algorithm **/
double **jacobi_algo() {
    int max_iter = 100, count = 0, max_i, max_j, check, i, j, m;
    double **P, **V, **A, **new_A, **result, eps = 1.0e-15, cal_eps, c, s;
    V = generate_matrix();
    P = generate_matrix();
    new_A = generate_matrix();
    A = laplacian_Lnorm();
    set_identity(V);
    do {
        set_identity(P);
        check = set_data(A, &c, &s, &max_i, &max_j);
        if (check == -1) {
            break;
        }
        P[max_i][max_i] = P[max_j][max_j] = c;
        P[max_i][max_j] = s;
        P[max_j][max_i] = -s;
        multiply(V, P);

        cal_newA(new_A, A, c, s, max_i, max_j);
        cal_eps = get_of_f_A(A) - get_of_f_A(new_A);
        A = new_A;
        count++;
    } while (count < max_iter && cal_eps > eps);

    result = (double **) calloc(n + 1, sizeof(double *));
    for (i = 0; i < n; ++i) {
        result[i] = (double *) calloc(n, sizeof(double));
    }
    for (j = 0; j < n + 1; ++j) {
        for (m = 0; m < n; ++m) {
            if (j == 0) {
                result[j][m] = A[m][m];
            } else {
                result[j][m] = V[j - 1][m];
            }
        }
    }
    free_mat(P);
    free_mat(V);
    free_mat(A);
    return result;
}

/** compare between doubles **/
int compare(const void *obj1, const void *obj2) {
    if (obj1 > obj2) {
        return 1;
    }
    if (obj1 < obj2) {
        return -1;
    }
    return 0;
}

/** calculate k (number of clusters) out of the eigengap **/
int eigengap_heuristic() {
    int i, max_i = 0, k;
    double **jacobi_mat, *eigenVals, max_val = 0, delta;
    jacobi_mat = jacobi_algo();
    eigenVals = jacobi_mat[0];
    qsort(eigenVals, n, sizeof(double), compare);
    for (i = 0; i < (int) (n / 2); ++i) {
        delta = fabs(eigenVals[i] - eigenVals[i + 1]);
        if (delta > max_val) {
            max_val = delta;
            max_i = i;
        }
    }
    for (k = 0; k < n + 1; ++k) {
        free(jacobi_mat[k]);
    }
    free(jacobi_mat);
    return max_i;
}


/** print a given matrix
void print_matrix(double **matrix) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (j != n - 1) {
                printf("%f , ", matrix[i][j]);
            } else {
                printf("%f ", matrix[i][j]);
            }
        }
        printf("\n");
    }
}
**/
