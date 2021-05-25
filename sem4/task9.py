import numpy as np
import matplotlib.pyplot as plt

class Task9:
    def eiler_method(func_, _range_, init_v, accuracy):
        n = 2
        h = (_range_[1] - _range_[0]) / n
        points = [[init_v[0]], [init_v[1]]]
        Task9.eiler_points_append(func_, h, n, points)
        diff = accuracy * 4
        while accuracy < diff:
            points_temp = points.copy()
            n *= 2
            h /= 2
            points = [[init_v[0]], [init_v[1]]]
            Task9.eiler_points_append(func_, h, n, points)
            diff = Task9.derivative(func_, h, points_temp[1], points_temp[0])
        return (n, diff, points)

    def eiler_method_mod(func_, _range_, init_v, accuracy):
        n = 2
        h = (_range_[1] - _range_[0]) / n
        points = [[init_v[0]], [init_v[1]]]
        Task9.modified_euler_points_append(func_, h, n, points)
        diff = accuracy * 4
        while diff > accuracy:
            points_temp = points.copy()
            n *= 2
            h /= 2
            points = [[init_v[0]], [init_v[1]]]
            Task9.modified_euler_points_append(func_, h, n, points)
            diff = Task9.derivative(func_, h, points_temp[1], points_temp[0])
        return (n, diff, points)

    def eiler_points_append(func_, h, n, points):
        for i in range(int(n)):
            points[0].append(points[0][i] + h)
            points[1].append(points[1][i] + h * func_(points[0][i], points[1][i]))

    def runge_method(func_, _range_, init_v, accuracy):
        n = 2
        h = (_range_[1] - _range_[0]) / n
        points = [[init_v[0]], [init_v[1]]]
        Task9.runge_points_append(func_, h, n, points)
        diff = accuracy * 4
        while diff > accuracy:
            points_temp = points.copy()
            n *= 2
            h /= 2
            points = [[init_v[0]], [init_v[1]]]
            Task9.runge_points_append(func_, h, n, points)
            diff = Task9.derivative(func_, h, points_temp[1], points_temp[0])
        return (n, diff, points)

    def runge_points_append(func_, h, n, points):
        for i in range(int(n)):
            points[0].append(points[0][i] + h)
            tmp = []
            tmp.append(h * func_(points[0][i], points[1][i]))
            tmp.append(h * func_(points[0][i] + h / 2, points[1][i] + tmp[0] / 2))
            tmp.append(h * func_(points[0][i] + h / 2, points[1][i] + tmp[1] / 2))
            tmp.append(h * func_(points[0][i] + h, points[1][i] + tmp[2]))
            points[1].append(points[1][i] + (1 / 6) * (tmp[0] + tmp[3] + 2 * (tmp[1] + tmp[2])))


    def modified_euler_points_append(func_, h, n, points):
        for i in range(int(n)):
            points[0].append(points[0][i] + h)
            points[1].append(points[1][i] +
                             h * func_(points[0][i] +
                                       h / 2, points[1][i] +
                                       (h / 2) * func_(points[0][i], points[1][i])))

    def derivative(f, h, old_data_y, old_data_x):
        size = len(old_data_y)
        d_y = abs(np.array([old_data_y[i] - (old_data_y[i] + h * f(old_data_x[i], old_data_y[i])) for i in range(1, size)]))
        return d_y.max()
def test1():
    f = lambda x, y: x**3 + 3*y
    initial_value = (0, 0)
    tol = 0.001
    res = Task9.eiler_method(f, [0, 1], initial_value, tol)
    print("\nEuler method: " + str(res[1]))
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.eiler_method_mod(f, [0, 1], initial_value, tol)
    print("\nModified Euler method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.runge_method(f, [0, 1], initial_value, tol)
    print("\nRunge method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

def test2():
    f = lambda x, y: np.sin(x) * np.cos(y)
    initial_value = (0, 0)
    tol = 0.001
    res = Task9.eiler_method(f, [0, 1], initial_value, tol)
    print("\nEuler method: " + str(res[1]))
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.eiler_method_mod(f, [0, 1], initial_value, tol)
    print("\nModified Euler method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.runge_method(f, [0, 1], initial_value, tol)
    print("\nRunge method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()
def test3():
    f = lambda x, y: x - y
    initial_value = (0, 0)
    tol = 0.001
    res = Task9.eiler_method(f, [0, 1], initial_value, tol)
    print("\nEuler method: " + str(res[1]))
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.eiler_method_mod(f, [0, 1], initial_value, tol)
    print("\nModified Euler method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()

    res = Task9.runge_method(f, [0, 1], initial_value, tol)
    print("\nRunge method: ")
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    plt.show()
if __name__ == "__main__":
    # f = lambda x, y: (0.9 * (1 - y ** 2)) / ((1 + 1.5) * x ** 2 + y ** 2 + 1)
    # initial_value = (0, 0)
    # tol = 0.001
    # res = Task9.eiler_method(f, [0, 1], initial_value, tol)
    # print("\nEuler method: " + str(res[1]))
    # print('Segments: ' + str(res[0]))
    # print('Error: ' + str(res[1]))
    # plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    # plt.show()
    #
    # res = Task9.eiler_method_mod(f, [0, 1], initial_value, tol)
    # print("\nModified Euler method: ")
    # print('Segments: ' + str(res[0]))
    # print('Error: ' + str(res[1]))
    # plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    # plt.show()
    #
    # res = Task9.runge_method(f, [0, 1], initial_value, tol)
    # print("\nRunge method: ")
    # print('Segments: ' + str(res[0]))
    # print('Error: ' + str(res[1]))
    # plt.scatter(res[2][0], res[2][1], s=0.5, color='b')
    #
    # plt.show()
    test3()
