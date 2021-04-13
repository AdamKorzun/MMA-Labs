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


if __name__ == '__main__':
    x_array = np.arange(0, 1.1, 0.1)
    y_array = Task6.get_array(x_array, 8, 4)


    approximation = np.poly1d(np.polyfit(x_array, y_array, 10))
    print(f"Интерполяционный многочлен Лагранжа: \n{Task6.lagrange_method(x_array, y_array)}")
    print(f"Интерполяционный многочлен Ньютона: \n{Task6.newton_method(x_array, y_array)}")
    print(f"Многочлен наилучшего приближения: \n{approximation}")
    print(f"Значение интерполяционного многочлена в 0.47: {Task6.lagrange_method(x_array, y_array)(0.47)}")
    print(f"Значение многочлена наилучшего приближения в 0.47: {approximation(0.47)}")
    y_plot = Task6.lagrange_method(x_array, y_array)(x_array)
    plt.plot(x_array, y_array, 'o', x_array, y_plot)
    plt.grid(True)
    plt.show()
    '''
    print("\nТестовый пример 1")
    x = np.array([-1.5, -0.75, 0, 0.75, 1.5])
    y = np.array([-14.1014, -0.931596, 0, 0.931596, 14.1014])
    test_lagrange = Task6lagrange(x, y)
    test_newton = newton(x, y)
    approximation = np.poly1d(np.polyfit(x, y, 3))
    print(f"Интерполяционный многочлен Лагранжа: \n{test_lagrange}")
    print(f"Интерполяционный многочлен Ньютона: \n{test_newton}")
    print(f"Многочлен наилучшего приближения: \n{approximation}")
    print(f"Значение интерполяционного многочлена в 0.47: {test_lagrange(0.47)}")
    print(f"Значение многочлена наилучшего приближения в 0.47: {approximation(0.47)}")
    x_new = np.arange(-2, 2, 0.1)
    y_plot = test_newton(x_new)
    plt.plot(x, y, 'o', x_new, y_plot)
    plt.grid(True)
    plt.show()

    print("\nТестовый пример 2")
    x = np.array([-5.0, -1.0, 0.0, 2.0])
    y = np.array([-2.0, 6.0, 1.0, 3.0])
    test_lagrange2 = lagrange(x, y)
    test_newton2 = newton(x, y)
    approximation = np.poly1d(np.polyfit(x, y, 3))
    print(f"Интерполяционный многочлен Лагранжа: \n{test_lagrange2}")
    print(f"Интерполяционный многочлен Ньютона: \n{test_newton2}")
    print(f"Многочлен наилучшего приближения: \n{approximation}")
    print(f"Значение интерполяционного многочлена в 0.47: {test_lagrange2(0.47)}")
    print(f"Значение многочлена наилучшего приближения в 0.47: {approximation(0.47)}")
    x_new = np.arange(-5.1, 2.1, 0.1)
    y_plot = test_newton2(x_new)
    plt.plot(x, y, 'o', x_new, y_plot)
    plt.grid(True)
    plt.show()
    '''
