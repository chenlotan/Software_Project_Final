#include "spkmeans.h"

int check_allocation_array(const double *p);
int check_allocation_2d_array(const double **p);


void print_matrix(double **matrix, int rows, int cols){
    int i,j;
    for(i=0; i<rows; ++i){
        for (int j = 0; j < cols; ++j) {
            if (j != cols - 1) {
                printf("%0.4f,", matrix[i][j]);
            } else {
                printf("%0.4f ", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

/** free memory of nxn matrix **/
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
        check_allocation_array(matrix[i]);
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
    printf("Weighted Matrix\n");
    print_matrix(w_mat, n, n);
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
    printf("ddg Matrix\n");
    print_matrix(diag_mat, n, n);
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
    printf("Lnorm Matrix\n");
    print_matrix(Lnorm_mat, n, n);
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
/** get symmetric matrix and return matrix of eigenvectors, when the first row is the eigenvalues **/
double **jacobi_algo(double **A, int n) {
    int max_iter = 100, count = 0, max_i, max_j, check, i, j, m;
    double **P, **V, **new_A, **result, eps = 1.0e-15, cal_eps, c, s;
    V = generate_matrix(n, n);
    P = generate_matrix(n, n);
    new_A = generate_matrix(n, n);
//    A = laplacian_Lnorm(data_points, n, dimension);
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
        check_allocation_array(result[i]);
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
    free_mat(P, n);
    free_mat(V, n);
    free_mat(A, n);
    printf("Jacobi Matrix\n");
    print_matrix(result, n,n);
    return result;
}


/** compare between doubles **/
int compare(const void *p1, const void *p2) {
    double obj1 = *(double *)p1;
    double obj2 = *(double *)p2;
    if (obj1 > obj2) {
        return 1;
    }
    if (obj1 < obj2) {
        return -1;
    }
    return 0;
}

/** calculate k (number of clusters) out of the eigengap **/
int eigengap_heuristic(double **data_points, int n, int dimension) {
    int i, max_i = 0, k;
    double **jacobi_mat, *eigenVals, max_val = 0, delta;
    jacobi_mat = jacobi_algo(data_points, n);
    eigenVals = jacobi_mat[0];

    qsort(eigenVals, n, sizeof(double), compare);
    for (k=0; k<n; k++){
        printf("%f, ", eigenVals[k]);
    }
    printf("\n");
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
/**transpose matrix**/
double **transpose(double **matrix, int rows, int cols){
    int i, j;
    double **result = (double **) malloc(cols * sizeof(double *));
    check_allocation_2d_array(result);
    for(i=0; i<cols; i++){
        result[i] = (double *) calloc(sizeof (double ), rows);
        check_allocation_array(result[i]);
        for (j=0; j<rows; j++){
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

/**comparator for qsort - in order to sort vectors**/
int compare_vec(void *v1, void *v2){
    double *vec1 = (double *)v1;
    double *vec2 = (double *)v2;
    if (vec1[0]>vec2[0]){
        return 1;
    }
    if (vec1[0]<vec2[0]){
        return -1;
    }
    return 0;
}

/**create matrix from eigenvectors of k minimal eigenvalues**/
double **k_min_eigenvectors(double **matrix, int k, int n){
    int i, j;
    double **transpose_mat, **result = (double **) malloc(sizeof (double *)*n);
    check_allocation_2d_array(transpose_mat);
    transpose_mat = transpose(matrix, n+1, n);
    qsort(transpose_mat, n, sizeof(double *), compare_vec);
    for(i=0; i<n; i++){
        result[i] = (double *) calloc(k, sizeof(double));
        check_allocation_array(result[i]);
        for(j=0; j<k; j++){
            result[i][j] = transpose_mat[j][i+1];
        }
    }
    free_mat(transpose_mat, n);
    return result;
}



/**create T matrix - vectors of k min eigenvalues normalized
 * in case of row of zeros - the row will stay the same**/
double **create_T_matrix(double **matrix, int k, int n){
    int i, j, sum_row;
    double  **result = k_min_eigenvectors(matrix, k, n);
    for(i=0; i<n; i++){
        sum_row=0;
        for(j=0; j<k; j++){
            sum_row += pow(matrix[i][j], 2);
        }
        for(j=0; j<k; j++){
            if (sum_row != 0){
                result[i][j] = result[i][j]/(sqrt(sum_row));
            }
        }
    }
    return result;
}



/** run full spkmeans algorithm. return T matrix, the rest of the algorithm will run from Python**/
double **sp_kmeans(double **data_points, int n, int dimension, int k) {
    double **T_matrix, **eigen_mat, **lnorm_mat;
    lnorm_mat = laplacian_Lnorm(data_points, n, dimension);
    eigen_mat = jacobi_algo(lnorm_mat, n);
    if (k==0){
        k = eigengap_heuristic(data_points, n, dimension);
    }

    T_matrix = create_T_matrix(eigen_mat, k, n);
    free_mat(eigen_mat, n+1);
    return T_matrix;
}

int check_allocation_array(const double *p) {
    if (p == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    return 0;
}

int check_allocation_2d_array(const double **p){
    if (p == NULL) {
        printf("An Error Has Occurred");
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
    check_allocation_array(res);
    ch = fscanf(file, "%s", buff);
    while((ch != 'n') && (ch != EOF)){
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



/***********************************handles only txt need to add support for csv files**********************************************/
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
            check_allocation_array(vector);
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
        printf("Invalid Input!");
        exit(1);
    }

}


int main(int argc, char *argv[]) {
    int goal, dimension, n, i, j;
    int *shape;
    char *file_name;
    double **data_points, **result;
    if (argc != 3){
        printf("Invalid Input!");
        exit(1);
    }
    goal = (int) strtol(argv[1], NULL, 10);
    file_name = argv[2];
    shape = find_shape(file_name);
    n = shape[0];
    dimension = shape[1];
    free(shape);
    data_points = read_file(file_name, n, dimension);
    result = sp_kmeans(data_points, n, dimension, 4);
    printf("Result\n");
    print_matrix(result, n, 4);
    goal = 10;
    switch (goal) {
        case 1:
            result = weight_matrix(data_points, n, dimension);
            print_matrix(result, n, n);
            break;
        case 2:
            result = diagonal_d_matrix(data_points, n, dimension);
            print_matrix(result, n,n);
            break;
        case 3:
            result = laplacian_Lnorm(data_points, n, dimension);
            print_matrix(result, n, n);
            break;
        case 4:
            result = jacobi_algo(data_points, n);
            print_matrix(result, n+1, n);
            break;
        default:
            printf("Invalid_Input!");
    }
    return 0;
}
