import numpy as np


class Task2:
    def get_matrix(variant):
        c_matrix = np.array([[0.01, 0, -0.02, 0, 0],
                             [0.01, 0.01, 0.02, 0, 0],
                             [0, 0.01, 0.01, 0, -0.02],
                             [0, 0, 0.01, 0.01, 0],
                             [0, 0, 0, 0.01, 0.01]])
        d_matrix = np.array([[1.33, 0.21, 0.17, 0.12, -0.13],
                             [-0.13, -1.33, 0.11, 0.17, 0.12],
                             [0.12, -0.13, -1.33, 0.11, 0.17],
                             [0.17, 0.12, -0.13, -1.33, 0.11],
                             [0.11, 0.67, 0.12, -0.13, -1.33]])
        #b_vector = np.array([[1.2], [2.2], [4.0], [0.0], [-1.2]])
        #res_matrix = np.concatenate((a_matrix, b_vector), axis=1)
        #return res_matrix
        a_matrix = c_matrix * variant + d_matrix
        return a_matrix

    def iteration_method():
        None

    def seidel(A, b, x, N, tol):
        maxIterations = 1000000
        xprev = [0.0 for i in range(N)]
        for i in range(maxIterations):
            for j in range(N):
                xprev[j] = x[j]
            for j in range(N):
                summ = 0.0
                for k in range(N):
                    if (k != j):
                        summ = summ + A[j][k] * x[k]
                x[j] = (b[j] - summ) / A[j][j]
            diff1norm = 0.0
            oldnorm = 0.0
            for j in range(N):
                diff1norm = diff1norm + abs(x[j] - xprev[j])
                oldnorm = oldnorm + abs(xprev[j])
            if oldnorm == 0.0:
                oldnorm = 1.0
            norm = diff1norm / oldnorm
            if (norm < tol) and i != 0:
                x_list = []
                for j in range(N):
                    x_list.append(x[j][0])
                return x_list



if __name__ == "__main__":
    a_matrix = Task2.get_matrix(8)
    b_vector = np.array([[1.2], [2.2], [4.0], [0.0], [-1.2]])
    guess = [0.0, 0.0,0.0, 0.0, 0.0]
    print(Task2.seidel(a_matrix.copy(), b_vector.copy(), guess, len(guess), 0.0001))
    #print(Task2.get_matrix(8))
