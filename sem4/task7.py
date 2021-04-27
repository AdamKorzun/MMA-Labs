import matplotlib.pyplot as plt
import numpy as np


class Task7:
    def interpolate(inp, point):
        nodes = np.array(inp[0])
        f = np.array(inp[1])

        if not np.any(np.diff(nodes) < 0):
            pass
        else:
            idx = np.argsort(nodes)
            nodes = nodes[idx]
            f = f[idx]

        x_diff = np.diff(nodes)
        y_diff = np.diff(f)

        n = len(nodes)

        arr_b = np.zeros(shape=(n, 1))
        arr = np.zeros(shape=(n, n))

        arr[0, 0] = 1
        arr[n - 1][n - 1] = 1
        lin_solution = Task7.lin(arr, arr_b, x_diff, y_diff, nodes)
        arr_d = np.zeros(shape=(n - 1, 1))
        arr_b = np.zeros(shape=(n - 1, 1))



        for i in range(0, len(arr_d)):
            arr_d[i] = (lin_solution[i + 1] - lin_solution[i]) / (
            3 * x_diff[i])
            arr_b[i] = (y_diff[i] / x_diff[i]) - (x_diff[i] / 3) * (
            2 * lin_solution[i] + lin_solution[i + 1])

        res = 0
        res = Task7.cidx(nodes, point, res)
        return f[res] + arr_b[res] * (
        point - nodes[res]) + lin_solution[res] * (
        point - nodes[res]) ** 2 + arr_d[res] * (
                point - nodes[res]) ** 3

    def cidx(nodes, point, res):
        for i in range(len(nodes) - 1):

            if nodes[i] <= point <= nodes[i + 1]:
                res = i
                break
        return res

    def lin(inp_array, b, dx, dy, nodes):

        for i in range(1, len(nodes) - 1):
            inp_array[i, i - 1] = dx[i - 1]
            inp_array[i, i + 1] = dx[i]
            inp_array[i, i] = 2 * (dx[i - 1] + dx[i])
            b[i, 0] = 3 * (dy[i] / dx[i] - dy[i - 1] / dx[i - 1])
        c = np.linalg.solve(inp_array, b)
        return c

def test_1():
    f = lambda x: x**2
    func = [[(x) for x in range(10)], [f(x)  for x in range(10)]]
    print(func)
    print(f(2.5))
    print(Task7.interpolate(func, 2.5))

def test_2():
    f = lambda x: x ** 2 * np.sin(x)
    func = [[(x) for x in range(10)], [f(x)  for x in range(10)]]
    print(func)
    print(f(2.5))
    print(Task7.interpolate(func, 2.5))

def test_3():
    f = lambda x: np.cos(x) * np.sin(x)
    func = [[(x) for x in range(1, 10)], [f(x)  for x in range(1, 10)]]
    print(func)
    print(f(2.5))
    print(Task7.interpolate(func, 2.5))

if __name__ == "__main__":
    f = lambda x : np.tanh(x)
    func = [[2 / 6 * x for x in range(6)], [f(2 / 6 * x) for x in range(6)]]
    print(Task7.interpolate(func, 1))
    #test_3()
    plt.scatter(func[0], func[1])
    plt.show()
