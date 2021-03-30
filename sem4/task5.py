import numpy as np
import math
import random

class Task5:
    def get_u_matrix(A):
        n = np.shape(A)[0]
        U = np.eye(n)
        U[(0,0)] = np.sqrt(A[(0,0)])
        U[0,1:] = [A[0][j] / U[0][0] for j in range(1, n)]

        for i in range(n):
            ss = sum([U[k][i]**2 for k in range(i)])
            U[i][i] = np.sqrt(A[i][i] - ss)
            for j in range(i+1, n):
                ss = sum([U[k][i] * U[k][j] for k in range(i)])
                U[i][j] = (A[i][j] - ss) / U[i][i]
        return U


    def eigen(a_matrix: np.array, tol):
        a_length = len(a_matrix[0])
        #A = A.copy()
        V = np.eye(a_length)
        rad = 0
        Uv = np.eye(a_length)
        valid = True
        while valid:
            i, j = 0, 1
            valid = False
            for ii in range(a_length):
                for jj in range(ii+1, a_length):
                    if abs(a_matrix[i][j]) <= abs(a_matrix[ii][jj]) and abs(np.real(a_matrix[ii][jj])) > tol:
                        i, j = ii, jj
                        valid = True

            U = np.eye(a_length)
            if a_matrix[j][j]-a_matrix[i][i] == 0:
                rad = np.pi/ 4
            else:
                rad = 0.5 * math.atan(2 * a_matrix[i][j] / (a_matrix[j][j]-a_matrix[i][i]))
            s = np.sin(rad)
            c = np.cos(rad)
            U[i][i] = U[j][j] = c
            U[i][j] = s
            U[j][i] = -s
            K = np.dot(np.transpose(U), a_matrix)
            a_matrix = np.dot(K, U)
            Uv = np.dot(Uv,U)
        eiges = [round(np.real(a_matrix[i][i]), 8) for i in range(a_length)]
        return (eiges, Uv)


def get_a_matrix(k):
    c_matrix = np.array([[0.2, 0.0, 0.2, 0.0, 0.0],
                         [0.0, 0.2, 0.0, 0.2, 0.0],
                         [0.2, 0.0, 0.2, 0.0, 0.2],
                         [0.0, 0.2, 0.0, 0.2, 0.0],
                         [0.0, 0.0, 0.2, 0.0, 0.2]])
    d_matrix = np.array([[2.33, 0.81, 0.67, 0.92, -0.53],
                         [0.81, 2.33, 0.81, 0.67, 0.92],
                         [0.67, 0.81, 2.33, 0.81, 0.92],
                         [0.92, 0.67, 0.81, 2.33, -0.53],
                         [-0.53, 0.92, 0.92, -0.53, 2.33]])
    return k * c_matrix + d_matrix
def test_eigen(matrix, tol):
    eigen_value, eigen_vector = Task5.eigen(matrix.copy(), tol)
    print(eigen_value)
    print(eigen_vector)
    verify = np.linalg.eig(matrix)
    print("\nEigen values (numpy):\n", verify[0])
    print("Eigen vectors (numpy):\n", verify[1])

if __name__ == "__main__":
    tol = 0.0001
    matrix = np.array([[-1, -6],
                       [2, 6]])
    test_eigen(matrix, tol)
    matrix = np.array([[5, 6, 3],
                       [-1, 0, 1],
                       [1, 2, -1]])
    test_eigen(matrix, tol)

    a = get_a_matrix(8)
    U = Task5.get_u_matrix(a)
    a_ = U.T @ U
    tol = 0.0001
    eiges, V = Task5.eigen(a_.copy(), tol)

    print("Eigen valuse:\n", eiges)
    print("Eigen vectors:\n", V)
    verify = np.linalg.eig(a)
    print("\nEigen values (numpy):\n", verify[0])
    print("Eigen vectors (numpy):\n", verify[1])
