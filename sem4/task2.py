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
        a_matrix = c_matrix * variant + d_matrix
        return a_matrix

    def simple_iterations(matrix, vector, eps):
        n = len(matrix)
        res_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    res_matrix[i][j] = -matrix[i][j] / matrix[i][i]
        res_vector = vector.copy()
        for i in range(n):
            res_vector[i] /= matrix[i][i]
        x = res_vector.copy()
        next_x = x.copy()
        while True:
            x = next_x.copy()
            for i in range(n):
                next_x[i] = 0
                for j in range(n):
                    next_x[i] += x[j]*res_matrix[i][j]
                next_x[i] += res_vector[i]
            norm = abs(next_x - x)
            if norm.max() < eps:
                break
        return x

    def norm(matrix):
        norm = np.zeros(len(matrix))
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                norm[i] += abs(matrix[i][j])
        return norm.max()

    def seidel(A, b, x, N, tol):
        max_iterations = 1000000
        xprev = [0.0 for i in range(N)]
        for i in range(max_iterations):
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
    print(a_matrix)
    np.set_printoptions(16)
    b_vector = np.array([[1.2], [2.2], [4.0], [0.0], [-1.2]])
    guess = [0.0, 0.0,0.0, 0.0, 0.0]
    tol = 0.00001
    print(Task2.seidel(a_matrix.copy(), b_vector.copy(), guess, len(guess), tol))
    print(Task2.simple_iterations(a_matrix.copy(), b_vector.copy(), tol))
    '''
    test_1 = np.array([[1, -1],
                       [0, 1]])
    b = np.array([[-5], [1]])
    guess =[0.0, 0.0]
    print(Task2.seidel(test_1, b, guess, len(guess) , 0.0001))
    '''
