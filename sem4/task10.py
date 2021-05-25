import numpy as np
import matplotlib.pyplot as plt

class Task10:
    def add_points(f, initial_value, h):
        points = [[initial_value[0]], [initial_value[1]]]
        points[0].append(points[0][0] + h)
        k = []
        k.append(h * f(points[0][0], points[1][0]))
        k.append(h * f(points[0][0] + h / 2, points[1][0] + k[0] / 2))
        k.append(h * f(points[0][0] + h / 2, points[1][0] + k[1] / 2))
        k.append(h * f(points[0][0] + h, points[1][0] + k[2]))
        points[1].append(points[1][0] + (1 / 6) * (k[0] + k[3] + 2 * (k[1] + k[2])))
        return points

    def get_diff(f, h, y_old, x_old):
        size = len(y_old)
        dY = np.array([y_old[i] - (y_old[i] + h * f(x_old[i], y_old[i])) for i in range(1, size)])
        dY = abs(dY)
        return dY.max()

    def adams(f, _range_, initial_value, tol):
        n = 2
        h = (_range_[1] - _range_[0]) / n
        points = Task10.add_points(f, initial_value, h)
        for i in range(int(n) - 1):
            points[0].append(points[0][i + 1] + h)
            points[1].append(points[1][i + 1] + h * ((3 / 2) * f(points[0][i + 1], points[1][i + 1]) - (1 / 2) * f(points[0][i], points[1][i])))
        diff = tol * 4
        while diff > tol:
            points_temp = points.copy()
            n *= 2
            h /= 2
            points = Task10.add_points(f, initial_value, h)
            for i in range(int(n) - 1):
                points[0].append(points[0][i + 1] + h)
                points[1].append(points[1][i + 1] + h * ((3 / 2) * f(points[0][i + 1], points[1][i + 1]) - (1 / 2) * f(points[0][i], points[1][i])))
            diff = Task10.get_diff(f, h, points_temp[1], points_temp[0])
        return (n, diff, points)

def test1():
    f = lambda x,y : x**3 + 3*y
    initial_value = (0, 0)
    solution_range = [0, 1]
    res = Task10.adams(f, solution_range, initial_value, 1e-3)
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='r')
    plt.show()

def test2():
    f = lambda x,y : np.sin(x) * np.cos(y)
    initial_value = (0, 0)
    solution_range = [0, 1]
    res = Task10.adams(f, solution_range, initial_value, 1e-3)
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='r')
    plt.show()

def test3():
    f = lambda x,y : x - y
    initial_value = (0, 0)
    solution_range = [0, 1]
    res = Task10.adams(f, solution_range, initial_value, 1e-3)
    print('Segments: ' + str(res[0]))
    print('Error: ' + str(res[1]))
    plt.scatter(res[2][0], res[2][1], s=0.5, color='r')
    plt.show()
if __name__ == "__main__":
    # m = 1.5
    # a = 0.9
    # f = lambda x,y : (a * (1 - y ** 2)) / ((1 + m) * x ** 2 + y ** 2 + 1)
    # initial_value = (0, 0)
    # solution_range = [0, 1]
    # res = Task10.adams(f, solution_range, initial_value, 1e-3)
    # print('Segments: ' + str(res[0]))
    # print('Error: ' + str(res[1]))
    # plt.scatter(res[2][0], res[2][1], s=0.5, color='r')
    # plt.show()
    test3()
