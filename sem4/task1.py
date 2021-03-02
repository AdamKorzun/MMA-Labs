import numpy as np


class Task1:
    def gaussian_elimination_1(a_matrix):
        n = len(a_matrix)
        x = np.zeros(n)
        for i in range(n):
            if a_matrix[i][i] == 0.0:
                print(a_matrix)
                # swap row with non-zero one
                swapped = False
                for j in range(i+1, len(a_matrix)):
                    if a_matrix[j][i] != 0:
                        a_matrix[[j,i]] = a_matrix[[i, j]]
                        swapped = True
                        break

                if not swapped:

                    raise ValueError()

            for j in range(i+1, n):
                ratio = a_matrix[j][i]/a_matrix[i][i]
                for k in range(n+1):
                    a_matrix[j][k] = a_matrix[j][k] - ratio * a_matrix[i][k]
        x[n-1] = a_matrix[n-1][n]/a_matrix[n-1][n-1]
        for i in range(n-2,-1,-1):
            x[i] = a_matrix[i][n]
            for j in range(i+1,n):
                x[i] = x[i] - a_matrix[i][j]*x[j]
            x[i] = x[i]/a_matrix[i][i]
        return x

    def gaussian_elimination_2(a_matrix):
        for i in range(len(a_matrix[0]) - 1):
            max_el = abs(a_matrix[i][i])
            max_row = i
            for j in range(i, len(a_matrix  )):
                if (abs(a_matrix[j,i]) > max_el):
                    max_el = abs(a_matrix[j,i])
                    max_pos = j
            a_matrix[[i, max_row]] = a_matrix[[max_row,i]]
        return Task1.gaussian_elimination_1(a_matrix)

    def gaussian_elimination_3(a_matrix):
        # don't do row swap if a[i][i] == 0 and a[i] == max_row (fixed ? )
        for k in range(len(a_matrix)):
            max_el = abs(a_matrix[0,0])
            max_index = k
            for i in range(k, len(a_matrix)):
                for j in range(len(a_matrix[i]) - 1):
                    if (abs(a_matrix[i,j]) > max_el and a_matrix[i][i] != 0):
                        max_el = abs(a_matrix[i,j])
                        max_index = i
            a_matrix[[k, max_index]] = a_matrix[[max_index,k]]

            for j in range(k+1, len(a_matrix)):
                ratio = a_matrix[j][k]/a_matrix[k][k]
                for f in range(len(a_matrix)+1):
                    a_matrix[j][f] = a_matrix[j][f] - ratio * a_matrix[k][f]
            #print('max row index: {}, max el: {}'.format(max_index, max_el))
            #print(a_matrix)
        n = len(a_matrix)
        x = np.zeros(n)
        x[n-1] = a_matrix[n-1][n]/a_matrix[n-1][n-1]
        for i in range(n-2,-1,-1):
            x[i] = a_matrix[i][n]
            for j in range(i+1,n):
                x[i] = x[i] - a_matrix[i][j]*x[j]
            x[i] = x[i]/a_matrix[i][i]
        return x

    def get_a_matrix(variant):
        c_matrix = np.array([[0.2, 0, 0.2, 0, 0],
                             [0, 0.2, 0, 0.2, 0],
                             [0.2, 0, 0.2, 0, 0.2],
                             [0, 0.2, 0, 0.2, 0],
                             [0, 0, 0.2, 0, 0.2]])
        d_matrix = np.array([[2.33, 0.81, 0.67, 0.92, -0.53],
                              [-0.53, 2.33, 0.81, 0.67, 0.92],
                              [0.92, -0.53, 2.33, 0.81, 0.67],
                              [0.67, 0.92, -0.53, 2.33, 0.81],
                              [0.81, 0.67, 0.92, -0.53, 2.33]])
        b_vector = np.array([[4.2],[4.2],[4.2],[4.2],[4.2]])
        k_constant = variant
        c_matrix *= k_constant
        a_matrix=c_matrix +  d_matrix
        res_matrix = np.concatenate((a_matrix, b_vector), axis=1)
        return res_matrix
    def test_eliminations(matrix):
        res_1 = Task1.gaussian_elimination_1(matrix.copy())
        res_2 = Task1.gaussian_elimination_2(matrix.copy())
        res_3 = Task1.gaussian_elimination_3(matrix.copy())
        print('gaussian_elimination_1: {}'.format(res_1))
        print('gaussian_elimination_2: {}'.format(res_2))
        print('gaussian_elimination_3: {}'.format(res_3))



if __name__ == '__main__':
    np.set_printoptions(precision=6)
    a_matrix = Task1.get_a_matrix(8)
    print(a_matrix)
    print(Task1.gaussian_elimination_1(a_matrix.copy()))
    print(Task1.gaussian_elimination_2(a_matrix.copy()))
    print(Task1.gaussian_elimination_3(a_matrix.copy()))
    test_1 = np.array([[0, 1, 1],
                       [1, -1, -5]])
    Task1.test_eliminations(test_1)
    test_2 = np.array([[3, 2, -5, -1],
                       [2, -1, 3, 13],
                       [1, 2, -1, 9]])
    Task1.test_eliminations(test_2)
    test_3 = Task1.get_a_matrix(4)
    Task1.test_eliminations(test_3)
