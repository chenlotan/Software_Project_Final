#include "spkmeans.h"

int check_allocation_double_array(double *p);

int check_allocation_2d_array(double **p);

int check_allocation_int_array(int *p);

void print_matrix(double **matrix, int rows, int cols){
    int i,j;
    for(i=0; i<rows; ++i){
        for (j = 0; j < cols; ++j) {
            if (j != cols - 1) {
                printf("%0.4f,", matrix[i][j]);
            } else {
                printf("%0.4f", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

/** free memory of matrix **/
void free_mat(double **matrix, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

/** returns pointer to a matrix that initialized to zeros **/
double **generate_matrix(int rows, int cols) {
    int i;
    double **matrix;
    matrix = (double **) calloc(rows, sizeof(double *));
    check_allocation_2d_array(matrix);
    for (i = 0; i < rows; ++i) {
        matrix[i] = (double *) calloc(cols, sizeof(double));
        check_allocation_double_array(matrix[i]);
    }
    return matrix;
}

/** sets the given matrix to an identity matrix **/
void set_identity(double **mat, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        mat[i][i] = 1;
        for (j = i + 1; j < n; ++j) {
            mat[i][j] = mat[j][i] = 0;
        }
    }
}

/** given two initialized matrices at size n, copying the data from the first matrix to the second one **/
void copy_matrix(double **matrix, double **copy_matrix, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            copy_matrix[i][j] = matrix[i][j];
        }
    }
}

/** multiply matrix V in matrix P in place **/
void multiply(double **V, double **P, int n) {
    int i, j, k;
    double **A;
    A = generate_matrix(n, n);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k) {
                A[i][j] += (V[i][k] * P[k][j]);
            }
        }
    }
    copy_matrix(A, V, n);
    free_mat(A, n);
}

/** computes the Euclidean norm between two vectors **/
double compute_weight(double vec1[], double vec2[], int dimension) {
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        sum += (double) pow(vec1[i] - vec2[i], 2);
    }
    return exp(-sqrt(sum) / 2);
}

/** calculate the weight matrix out of the data points **/
double **weight_matrix(double **data_points, int n, int dimension) {
    int i, j;
    double **w_mat = generate_matrix(n, n);
    for (i = 0; i < n; ++i) {
        w_mat[i][i] = 0;
        for (j = i + 1; j < n; ++j) {
            w_mat[i][j] = compute_weight(data_points[i], data_points[j], dimension);
            w_mat[j][i] = w_mat[i][j];
        }
    }
    return w_mat;
}

/** computes the diagonal matrix out of the weight matrix **/
double **diagonal_d_matrix(double **data_points, int n, int dimension) {
    double **diag_mat, **w_mat, sum;
    int i, j;
    diag_mat = generate_matrix(n, n);
    w_mat = weight_matrix(data_points, n, dimension);
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
double **laplacian_Lnorm(double **data_points, int n, int dimension) {
    int i, j;
    double **diag_mat, **w_mat, **Lnorm_mat;

    w_mat = weight_matrix(data_points, n, dimension);
    diag_mat = diagonal_d_matrix(data_points, n, dimension);
    Lnorm_mat = generate_matrix(n, n);

    for (i = 0; i < n; ++i) {
        Lnorm_mat[i][i] = 1 - (w_mat[i][i] / diag_mat[i][i]);
        for (j = i + 1; j < n; ++j) {
            Lnorm_mat[i][j] = - w_mat[i][j] / sqrt(diag_mat[i][i] * diag_mat[j][j]);
            Lnorm_mat[j][i] = Lnorm_mat[i][j];
        }
    }
    free_mat(diag_mat, n);
    free_mat(w_mat, n);
    return Lnorm_mat;
}

/** calculate the f(A)^2 of a given matrix **/
double get_of_f_A(double **A, int n) {
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
void cal_newA(double **new_A, double **A, double c, double s, int max_i, int max_j, int n) {
    int r;
    for (r = 0; r < n; ++r) {
        if (r != max_i && r != max_j) {
            new_A[r][max_i] = new_A[max_i][r] = c * A[r][max_i] - s * A[r][max_j];
            new_A[r][max_j] = new_A[max_j][r] =c * A[r][max_j] + s * A[r][max_i];
        }
    }
    new_A[max_i][max_i] = pow(c, 2) * A[max_i][max_i] + pow(s, 2) * A[max_j][max_j] - 2 * s * c * A[max_i][max_j];
    new_A[max_j][max_j] = pow(c, 2) * A[max_j][max_j] + pow(s, 2) * A[max_i][max_i] + 2 * s * c * A[max_i][max_j];
    new_A[max_i][max_j] = new_A[max_j][max_i] = 0;
}

/** calculate c and s from the Lnorm matrix and returns 0 if it was valid and -1 elsewhere **/
int set_data(double **A, double *c, double *s, int *max_i, int *max_j, int n) {
    int i, j;
    double max_Aij = 0, theta, t;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            if (fabs(A[i][j]) > max_Aij) {
                max_Aij = fabs(A[i][j]);
                *max_i = i;
                *max_j = j;
            }
        }
    }
    if (max_Aij == 0) {
        return -1;
    }
    theta = (A[*max_j][*max_j] - A[*max_i][*max_i]) / (2 * A[*max_i][*max_j]);
    t = ((theta >= 0) - (theta < 0)) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    *c = 1 / sqrt(pow(t, 2) + 1);
    *s = t * (*c);
    return 0;
}

/** Jacobi Algorithm **/
double **jacobi_algo(double **A, int n) {
    int max_iter = 100, count = 0, max_i, max_j, check, i, j, m;
    double **P, **V, **new_A, **result, eps = 1.0e-5, cal_eps, c, s;
    V = generate_matrix(n, n);
    P = generate_matrix(n, n);
    new_A = generate_matrix(n, n);
    copy_matrix(A, new_A, n);
    set_identity(V, n);
    do {
        set_identity(P, n);
        check = set_data(A, &c, &s, &max_i, &max_j, n);
        if (check == -1) {
            break;
        }
        P[max_i][max_i] = P[max_j][max_j] = c;
        P[max_i][max_j] = s;
        P[max_j][max_i] = -s;
        multiply(V, P, n);

        cal_newA(new_A, A, c, s, max_i, max_j, n);
        cal_eps = get_of_f_A(A, n) - get_of_f_A(new_A, n);
        copy_matrix(new_A, A, n);
        count++;
    } while (count < max_iter && cal_eps > eps);

    result = (double **) calloc(n + 1, sizeof(double *));
    check_allocation_2d_array(result);
    for (i = 0; i < n + 1; ++i) {
        result[i] = (double *) calloc(n, sizeof(double));
        check_allocation_double_array(result[i]);
    }
    for (j = 0; j < n + 1; ++j) {
        for (m = 0; m < n; ++m) {
            if (j == 0) {
                if (-0.0001 < A[m][m] && A[m][m] < 0){
                    result[j][m] = -A[m][m];
                }
                else {
                    result[j][m] = A[m][m];
                }
            } else {
                result[j][m] = V[j - 1][m];
            }
        }
    }
    free_mat(P, n);
    free_mat(V, n);
    return result;
}

/** calculate k (number of clusters) out of the eigengap **/
int eigengap_heuristic(double *eigen_vals, int n) {
    int i, max_i = 0;
    double max_val = 0, delta;
    for (i = 0; i < (int) (n / 2); ++i) {
        delta = fabs(eigen_vals[i] - eigen_vals[i + 1]);
        if (delta > max_val) {
            max_val = delta;
            max_i = i;
        }
    }
    return max_i+1;
}

/**transpose matrix**/
double **transpose(double **matrix, int rows, int cols){
    int i, j;
    double **result = (double **) malloc(cols * sizeof(double *));
    check_allocation_2d_array(result);
    for(i=0; i<cols; i++){
        result[i] = (double *) calloc(sizeof (double ), rows);
        check_allocation_double_array(result[i]);
        for (j=0; j<rows; j++){
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

/**comparator for qsort - in order to sort vectors**/
int compare_vec(const void *v1, const void *v2){
    double *vec1 = *(double * const *)v1;
    double *vec2 = *(double * const *)v2;
    if (vec1[0]>vec2[0]){
        return 1;
    }
    if (vec1[0]<vec2[0]){
        return -1;
    }
    return 0;
}

/**transpose matrix and sort - return value is vector of the sorted eigenvalues**/
double** transpose_and_sort(double **jacobi_mat, int n){
    double **transposed;
    transposed = transpose(jacobi_mat, n+1, n);
    qsort(transposed, n, sizeof(double *), compare_vec);
    return transposed;
}

double *get_eigenVals(double **matrix, int n){
    int i;
    double *eigen_vals = (double *) calloc(n, sizeof (double));
    for (i = 0; i<n; i++){
        eigen_vals[i] = matrix[i][0];
    }
    return eigen_vals;
}


/**create T matrix - vectors of k min eigenvalues normalized
 * in case of row of zeros - the row will stay the same**/
double **create_T_matrix(double **matrix, int k, int n){
    int i, j;
    double sum_row;
    double **result;
    result = (double **) malloc(n*sizeof(double *));
    for(i=0; i<n; i++){
        result[i] = (double *) calloc(k, sizeof(double ));
        sum_row=0;
        for(j=0; j<k; j++){
            sum_row += pow(matrix[i+1][j], 2);
        }
        for(j=0; j<k; j++){
            if (sum_row != 0){
                result[i][j] = matrix[i+1][j]/(sqrt(sum_row));
            }
        }
    }
    return result;
}

int check_allocation_double_array(double *p) {
    if (p == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return 0;
}

int check_allocation_int_array(int *p) {
    if (p == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return 0;
}

int check_allocation_2d_array(double **p){
    if (p == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return 0;
}

/**--------------------------------------functions for running code from c file--------------------------------------------**/

int find_dimension(char line[]) {
    int i = 0;
    char *ptr = strtok(line, ",");
    while (ptr != NULL) {
        ptr = strtok(NULL, ",");
        i++;
    }
    return i;
}

/**read input file and find n and dimension**/
int *find_shape(char *fileName){
    FILE *file = fopen(fileName, "r");
    char buff[1024];
    int ch, n = 0;
    int *res = (int *)calloc(2,sizeof(int));
    check_allocation_int_array(res);
    if (file){
        ch = fscanf(file, "%s", buff);
        while((ch != '\n') && (ch != EOF)){
            if (n == 0){
                res[1] = find_dimension(buff);
            }
            n++;
            ch = fscanf(file, "%s", buff);
        }
        fclose(file);
        res[0] = n;
        return res;
    }
    else {
        printf("Invalid Input!\n");
        exit(1);
    }
    return n;
}

/**read file and create vectors list**/
double **read_file(char fileName[], int n, int dimension) {
    FILE *file = fopen(fileName, "r");
    char buff[1024], *ptr;
    int ch, i = 0;
    double *vector, *place;
    double **all_vectors;
    if (file) {
        all_vectors = (double **) malloc(n * sizeof(double *));
        check_allocation_2d_array(all_vectors);
        ch = fscanf(file, "%s", buff);
        while ((ch != '\n') && (ch != EOF)) {
            vector = (double *) (malloc(dimension * sizeof(double)));
            check_allocation_double_array(vector);
            place = vector;
            ptr = strtok(buff, ",");
            while (ptr != NULL) {
                *(vector) = strtod(ptr, NULL);
                ptr = strtok(NULL, ",");
                vector++;
            }
            all_vectors[i] = place;
            i++;
            ch = fscanf(file, "%s", buff);
        }
        fclose(file);
        return all_vectors;
    } else {
        printf("Invalid Input!\n");
        exit(1);
    }

}



int main(int argc, char *argv[]) {
    int dimension, n;
    int *shape;
    char *file_name, *goal;
    double **data_points, **result;
    if (argc != 3){
        printf("Invalid Input!\n");
        exit(1);
    }
    goal = argv[1];
    file_name = argv[2];
    shape = find_shape(file_name);
    n = shape[0];
    dimension = shape[1];
    free(shape);
    data_points = read_file(file_name, n, dimension);
    if (strcmp(goal, "wam") == 0) {
        result = weight_matrix(data_points, n, dimension);
        print_matrix(result, n, n);
        free_mat(result, n);
        free_mat(data_points, n);
        return 0;
    }
    if (strcmp(goal, "ddg") == 0) {
        result = diagonal_d_matrix(data_points, n, dimension);
        print_matrix(result, n, n);
        free_mat(result, n);
        free_mat(data_points, n);
        return 0;
    }
    if (strcmp(goal, "lnorm") == 0) {
        result = laplacian_Lnorm(data_points, n, dimension);
        print_matrix(result, n, n);
        free_mat(result, n);
        free_mat(data_points, n);
        return 0;
    }
    if (strcmp(goal, "jacobi") == 0) {
        result = jacobi_algo(data_points, n);
        print_matrix(result, n + 1, n);
        free_mat(result, n + 1);
        free_mat(data_points, n);
        return 0;
    }
    else printf("Invalid_Input!\n");
    return 0;
}
