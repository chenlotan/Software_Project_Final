# Todo
# 1. choose function by cases according to the argument
# the module should support calling the functions:
#   1.full spkmeans,
#   2.calculate weighted matrix,
#   3.calculate diagonal degree matrix,
#   4.normalized graph laplacian,
#   5.eigenvalues and eigenvectors (jacoby).
#   6.heuristic - to calculate k
# 2. write function to spkmeans - that runs all levels and gives final output
# 3. write C-API module code
# notes:
# as written in the forum (https://moodle.tau.ac.il/mod/forum/discuss.php?d=64372):
#   epsilon value should be 0
#   max iter should be 300

import numpy as np
import pandas as pd
import sys
import mykmeanssp

def validate(condition):
    if not condition:
        print('Invalid Input!')
        exit(1)


def check_if_float(num):
    try:
        float(num)
    except ValueError:
        print('Invalid Input!')
        exit(1)


def read_file_to_df(file):
    if file.split(".")[-1] == "txt":
        return pd.read_csv(file, sep=",", header=None)
    else:
        return pd.read_csv(file)


def initialize_centroids(vectors, k):
    mu = np.zeros((k, len(vectors[0])))
    mu_index = np.zeros(k)
    np.random.seed(0)
    N = len(vectors)
    # print(vectors.shape[1])
    mu_index[0] = np.random.choice(N)
    mu[0] = vectors[(int(mu_index[0]))]
    for i in range(1, k):
        D = np.apply_along_axis(lambda vector: min(np.linalg.norm(vector-mu[:i], axis=1)**2), axis=1, arr=vectors)
        p = D/np.sum(D)
        mu_index[i] = np.random.choice(N, p=p)
        mu[i] = vectors[int(mu_index[i])]
    return mu_index.astype(int), mu.astype(float)


def fix_jacobi_eigenVals(matrix):
    for j in range(len(matrix[0])):
        if matrix[0][i] == 0:
            matrix[0][i] = 0.0000


def print_matrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if j < len(matrix[0]) - 1:
                print(str("{0:.4f}".format(matrix[i][j])) + ",", end="")
            else:
                print(str("{0:.4f}".format(matrix[i][j])))

goals = {"spk": 1, "wam": 2, "ddg": 3, "lnorm": 4, "jacobi": 5}
args = sys.argv[1:]
validate(len(args) == 3)
k, goal, input_filename = args[0], args[1], args[2]
validate(args[0].isdigit())
k = int(k)
goal = goals[goal]
epsilon = 0
max_iter = 300
vectors = read_file_to_df(input_filename)
n = vectors.shape[0]
dimension = vectors.shape[1]
# validate((k < n) & (k != 1))
result = mykmeanssp.spk_ext(k, goal, n, dimension, vectors.to_numpy().tolist())

if goal == 5:
    fix_jacobi_eigenVals(result)

if goal == 1:
    if k == 0:
        k = len(result[0])
        print("k =", k)
    centroids_index, centroids = initialize_centroids(result, k)
    final_centroids = mykmeanssp.fit(k, k, n, max_iter, epsilon, centroids.tolist(), result)
    for i in range(len(centroids_index)):
        if i < len(centroids_index) - 1:
            print(str(centroids_index[i]) + ",", end="")
        else:
            print(str(centroids_index[i]))
    print(final_centroids)
    print_matrix(final_centroids)

    # for i in range(len(final_centroids)):
    #     for j in range(len(final_centroids[0])):
    #         if j < len(final_centroids[0]) - 1:
    #             print(str("{0:.4f}".format(final_centroids[i][j])) + ",", end="")
    #         else:
    #             print(str("{0:.4f}".format(final_centroids[i][j])))


else:
    print_matrix(result)


