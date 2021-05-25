import numpy as np
import matplotlib.pyplot as plt


def eiler_method(func_, _range_, init_v, accuracy):
    n = 2
    h = (_range_[1] - _range_[0]) / n
    points = [[init_v[0]], [init_v[1]]]
    eiler_points_append(func_, h, n, points)
    diff = accuracy * 4
    while accuracy < diff:
        points_temp = points.copy()
        n *= 2
        h /= 2
        points = [[init_v[0]], [init_v[1]]]
        eiler_points_append(func_, h, n, points)
        diff = get_derivative(func_, h, points_temp[1], points_temp[0])
    print("\nEuler method: ")
    print("Number of segments: " + str(n))
    print("Error: " + str(diff))
    plt.scatter(points[0], points[1], s=0.5, color='blue')
    plt.show()


def eiler_modified_method(func_, _range_, init_v, accuracy):
    n = 2
    h = (_range_[1] - _range_[0]) / n
    points = [[init_v[0]], [init_v[1]]]
    modified_euler_points_append(func_, h, n, points)
    diff = accuracy * 4
    while diff > accuracy:
        points_temp = points.copy()
        n *= 2
        h /= 2
        points = [[init_v[0]], [init_v[1]]]
        modified_euler_points_append(func_, h, n, points)
        diff = get_derivative(func_, h, points_temp[1], points_temp[0])
    print("\nModified Euler method: ")
    print("Number of segments: " + str(n))
    print("Error: " + str(diff))
    plt.scatter(points[0], points[1], s=0.5, color='green')
    plt.show()


def eiler_points_append(func_, h, n, points):
    for i in range(int(n)):
        points[0].append(points[0][i] + h)
        points[1].append(points[1][i] + h * func_(points[0][i], points[1][i]))


def runge_method(func_, _range_, init_v, accuracy):
    n = 2
    h = (_range_[1] - _range_[0]) / n
    points = [[init_v[0]], [init_v[1]]]
    runge_points_append(func_, h, n, points)
    diff = accuracy * 4
    while diff > accuracy:
        points_temp = points.copy()
        n *= 2
        h /= 2
        points = [[init_v[0]], [init_v[1]]]
        runge_points_append(func_, h, n, points)
        diff = get_derivative(func_, h, points_temp[1], points_temp[0])
    print("\nRunge method: ")
    print("Number of segments: " + str(n))
    print("Error: " + str(diff))
    plt.scatter(points[0], points[1], s=0.5, color='red')
    plt.show()


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


target_f = lambda x, y: (1.3 * (1 - y ** 2)) / ((1 + 1) * x ** 2 + y ** 2 + 1)
initial_value = (0, 0)


def get_derivative(f, h, old_data_y, old_data_x):
    size = len(old_data_y)
    d_y = abs(np.array([old_data_y[i] - (old_data_y[i] + h * f(old_data_x[i], old_data_y[i])) for i in range(1, size)]))
    return d_y.max()


def main():
    eiler_method(target_f, [0, 1], initial_value, 1e-3)
    eiler_modified_method(target_f, [0, 1], initial_value, 1e-3)
    runge_method(target_f, [0, 1], initial_value, 1e-3)
    plt.show()


if __name__ == "__main__":
    main()
