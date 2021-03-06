
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
    mu_index[0] = np.random.choice(N)
    mu[0] = vectors[(int(mu_index[0]))]
    for i in range(1, k):
        D = np.apply_along_axis(lambda vector: min(np.linalg.norm(vector-mu[:i], axis=1)**2), axis=1, arr=vectors)
        p = D/np.sum(D)
        mu_index[i] = np.random.choice(N, p=p)
        mu[i] = vectors[int(mu_index[i])]
    return mu_index.astype(int), mu.astype(float)


def print_matrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if j < len(matrix[0]) - 1:
                print(str("{0:.4f}".format(matrix[i][j])) + ",", end="")
            else:
                print(str("{0:.4f}".format(matrix[i][j])))


EPSILON = 0
MAX_ITER = 300
goals = {"spk": 1, "wam": 2, "ddg": 3, "lnorm": 4, "jacobi": 5}
args = sys.argv[1:]
validate(len(args) == 3)
k, goal, input_filename = args[0], args[1], args[2]
validate(args[0].isdigit())
k = int(k)
goal = goals[goal]
vectors = read_file_to_df(input_filename)
n = vectors.shape[0]
dimension = vectors.shape[1]
validate((0 <= k < n) & (k != 1))
result = mykmeanssp.spk_ext(k, goal, n, dimension, vectors.to_numpy().tolist())

if goal == 1:
    if k == 0:
        k = len(result[0])
    centroids_index, centroids = initialize_centroids(result, k)
    final_centroids = mykmeanssp.fit(k, k, n, MAX_ITER, EPSILON, centroids.tolist(), result)
    for i in range(len(centroids_index)):
        if i < len(centroids_index) - 1:
            print(str(centroids_index[i]) + ",", end="")
        else:
            print(str(centroids_index[i]))
    print_matrix(final_centroids)

else:
    print_matrix(result)
