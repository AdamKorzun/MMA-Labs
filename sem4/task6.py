import numpy as np
import matplotlib.pyplot as plt

class Task6:

    def lagrange_method(x_array, f):
        result = np.poly1d(0)

        for j in range(len(f)):
            poly1 = np.poly1d(1)
            for i in range(len(x_array)):
                if i != j:
                    poly1 = np.polymul(poly1, np.poly1d([x_array[i]], True) / (
                    x_array[j] - x_array[i]))

            result = np.polyadd(result, f[j] * poly1)
        return result


    def newton_method(x_array, y_array):
        solution = np.poly1d(0)
        arr = np.copy(y_array)

        for i in range(1, len(x_array)):
            arr[i:len(x_array)] = (arr[i:len(x_array)] - arr[i - 1]) / (
            x_array[i:len(x_array)] - x_array[i - 1])

        for i in range(len(arr)):
            poly = np.poly1d(1)
            for j in range(i):
                poly *= np.poly1d([x_array[j]], True)
            poly *= arr[i]
            solution += poly
        return solution

    def get_array(x, k, m):
        p = np.array([0.0, 0.41, 0.79, 1.13, 1.46, 1.76, 2.04,
        2.3, 2.55, 2.79, 3.01])
        y = np.zeros(len(x))
        for i in range(len(x)):
            y[i] = p[i] + ((-1) ** k) * m
        return y

def test_1():
    x_array = np.array([1, 2, 3])
    y_array = np.array([1, 4, 9])
    print("Lagrange polynomial: \n{}".format(Task6.lagrange_method(
    x_array, y_array)))

    print("Newton polynomial: \n{}".format(Task6.newton_method(x_array, y_array)))
    print("Lagrange polynomial(1.5): {}".format(Task6.lagrange_method(
    x_array, y_array)(1.5)))
    print("Newton polynomial(1.5): {}".format(Task6.newton_method(x_array, y_array)(1.5)))
    plt.plot(x_array, y_array, 'o', x_array, Task6.lagrange_method(
    x_array, y_array)(x_array))
    plt.grid(True)
    plt.show()

def test_2():
    x_array = np.array([1, 2, 3, 4])
    y_array = np.array([1, 8, 27, 64])
    print("Lagrange polynomial: \n{}".format(Task6.lagrange_method(
    x_array, y_array)))

    print("Newton polynomial: \n{}".format(Task6.newton_method(x_array, y_array)))
    print("Lagrange polynomial(2.31): {}".format(Task6.lagrange_method(
    x_array, y_array)(2.31)))
    print("Newton polynomial(2.31): {}".format(Task6.newton_method(x_array, y_array)(2.31)))
    plt.plot(x_array, y_array, 'o', x_array, Task6.lagrange_method(
    x_array, y_array)(x_array))
    plt.grid(True)
    plt.show()

def test_3():
    x_array = np.array([-1, 1, 2])
    y_array = np.array([2, 3, -1])
    print("Lagrange polynomial: \n{}".format(Task6.lagrange_method(
    x_array, y_array)))

    print("Newton polynomial: \n{}".format(Task6.newton_method(x_array, y_array)))
    print("Lagrange polynomial(2.31): {}".format(Task6.lagrange_method(
    x_array, y_array)(2.31)))
    print("Newton polynomial(2.31): {}".format(Task6.newton_method(x_array, y_array)(2.31)))
    plt.plot(x_array, y_array, 'o', x_array, Task6.lagrange_method(
    x_array, y_array)(x_array))
    plt.grid(True)
    plt.show()

if __name__ == '__main__':

    x_array = np.arange(0, 1.1, 0.1)
    y_array = Task6.get_array(x_array, 8, 4)
    #x_array = np.array([i for i in range(0, 1.1, 0.1)])

    approximation = np.poly1d(np.polyfit(x_array, y_array, 10))
    print("Lagrange  polynomial: \n{}".format(Task6.lagrange_method(
    x_array, y_array)))
    print("Newton polynomial: \n{}".format(Task6.newton_method(x_array, y_array)))
    print("Approximation: \n{}".format(approximation))
    print("Lagrange polynomial(0.47): {}".format(Task6.lagrange_method(
    x_array, y_array)(0.47)))
    print("Newton polynomial(0.47): {}".format(Task6.newton_method(x_array, y_array)(0.47)))
    plt.plot(x_array, y_array, 'o', x_array, Task6.lagrange_method(
    x_array, y_array)(x_array))

    plt.grid(True)
    plt.show()
