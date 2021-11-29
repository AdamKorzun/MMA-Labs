import sympy
import numpy
import math
from matplotlib import pyplot
from prettytable import PrettyTable


x = sympy.Symbol('x')
start = -1
finish = 1
y_arrays = []
k = 0.2
T = 1
Y2 = []
g1 = lambda t: 0
g2 = lambda t: 0
fi = lambda x: x * ( 1 - x)
data_table = []

class Util:
    def solve(tau, h, t, last_layer = None):
        fx = 1
        A = []
        B = []
        Y = []
        i = 0
        if (True):
            pass
        xi = start
        while xi < finish - 1e-5:
            if t != 0:
                a = -k * tau
                b = h**2 + 2 * k * tau
                c = -k * tau
                d = tau * h**2 * fx + h**2 * last_layer[i]

                if i == 0:
                    Ai = - c / b
                    Bi = (d - a * g1(t)) / b
                else:
                    Ai = -c / (b + a * A[-1])
                    Bi = (d - a * B[-1]) / (b + a * A[-1])

                A.append(Ai)
                B.append(Bi)
            else:
                Y.append(fi(xi))
            xi += h
            i += 1

        if t == 0:
            Y.append(finish ** 2)
            return Y

        Y = [0] * 11

        Y[0] = g1(t)
        Y[-2] = (A[-2] * h * g2(t) + B[-2]) / (1 - A[-2])
        Y[-1] = Y[-2] + 4*h  # yb

        while i > 1:
            i -= 1
            Y[i] = A[i] * Y[i + 1] + B[i]

        return Y


    def solve2(tau, h, t, j):
        Y = [0] * 11
        fx = 1
        xi = start

        if j == 0:
            i = 0
            while xi < finish:
                Y[i] = fi(xi)
                i += 1
                xi += h
        else:
            Y[0] = g1(t)
            for i in range(1, 10):
                Y[i] = Y2[j - 1][i] + tau / h ** 2 * (Y2[j - 1][i + 1] - 2 * Y2[j - 1][i] + Y2[j - 1][i - 1]) #+ tau * fx
                xi += h

            Y[-2] = Y2[j - 1][-2] + tau / h ** 2 * (2 * h * g2(t) - Y2[j - 1][-2] + Y2[j - 1][-3]) #+ tau * fx
            Y[-1] = Y[-2] + 2 * h *g2(t)

        Y2.append(Y)


    def implicit_function(h, tau):
        y_arrays.clear()
        t = 0
        for i in range(0, 10):
            if i == 0:
                y_arrays.append(Util.solve(tau, h, t))
            else:
                y_arrays.append(Util.solve(tau, h, t, y_arrays[-1]))
            if t == 0:
                t = 2 * tau
            else:
                t *= 2
            if t > T:
                break


    def explicit_function(h, tau):
        t = 0
        Y2.clear()
        for i in range(0, 10):
            Util.solve2(tau, h, t, i)
            if t == 0:
                t = 2 * tau
            else:
                t *= 2
            if t > T:
                break


    def print_grafic(array, h, method, t):
        x = numpy.arange(start, finish + h, h)
        pyplot.plot(x, array, label=f't={t}')
        pyplot.legend()
        pyplot.title(method + " метод ")



    def create_table():
        table = PrettyTable(["h", "tau", "std(t=2tau)", "std(t=4tau)", "max(abs((t=2tau))", "max(abs(Mmod(t=4tau))"])
        for data in data_table:
            table.add_row(data)
        print(table)


    def my_round(value):
        return "{:.5f}".format(value)


def main():
    h = (finish - start) / 10
    tau = 0.5 * h ** 2 / k
    tau_test = tau
    N = h
    j = 1
    for i in range(0, 4):
        Util.implicit_function(h, tau_test)
        max_2_tau = 0
        for u in y_arrays[1]:
            if max_2_tau < u:
                max_2_tau = u

        max_4_tau = 0
        for u in y_arrays[2]:
            if max_4_tau < u:
                max_4_tau = u

        data_table.append([N, Util.my_round(tau_test), Util.my_round(y_arrays[1][j]), Util.my_round(y_arrays[2][j]), max_2_tau, max_4_tau])
        N *= 2
        j += 2
        tau_test /= 4

    Util.create_table()

    Util.implicit_function(h, tau)
    t = 0
    for array in y_arrays:
        Util.print_grafic(array, h, "Неявный", "{t:.9f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()

    data_table.clear()
    tau_test = tau
    N = h
    j = 1
    for i in range(0, 4):
        Util.explicit_function(h, tau_test)
        max_2_tau = 0
        for u in Y2[1]:
            if max_2_tau < u:
                max_2_tau = u

        max_4_tau = 0
        for u in Y2[2]:
            if max_4_tau < u:
                max_4_tau = u

        data_table.append([N, Util.my_round(tau_test), Util.my_round(Y2[1][j]), Util.my_round(Y2[2][j]), max_2_tau, max_4_tau])
        N *= 2
        j += 2
        tau_test /= 4

    Util.create_table()
    Util.explicit_function(h, tau)
    t = 0
    for array in Y2:
        Util.print_grafic(array, h, "Явный", "{t:.9f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()

    tau = h ** 2 / 6
    data_table.clear()
    tau_test = tau
    N = h
    j = 1
    for i in range(0, 4):
        Util.implicit_function(h, tau_test)
        max_2_tau = 0
        for u in y_arrays[1]:
            if max_2_tau < u:
                max_2_tau = u

        max_4_tau = 0
        for u in y_arrays[2]:
            if max_4_tau < u:
                max_4_tau = u

        data_table.append([N, Util.my_round(tau_test), Util.my_round(y_arrays[1][j]), Util.my_round(y_arrays[2][j]), max_2_tau, max_4_tau])
        N *= 2
        j += 2
        tau_test /= 4

    Util.create_table()
    Util.explicit_function(h, tau)
    t = 0
    for array in Y2:
        Util.print_grafic(array, h, "Явный (t=h^2 / 6)", "{t:.9f}".format(t=t))
        if t == 0:
            t = 2 * tau
        else:
            t *= 2
    pyplot.show()


if __name__ == "__main__":
    main()
