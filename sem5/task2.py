import sympy as sp
import numpy as np
from IPython.display import display
import math
from matplotlib import pylab

class Utils:
    def plot(functions, start_x, end_x, dx):
        for function in functions:
            dx = 0.01
            xlist = np.arange(start_x, end_x, dx)
            ylist = [function(p) for p in xlist]
            pylab.plot(xlist, ylist)
        pylab.grid(True)
        pylab.show()


    def get_basis_function(n):
        if not n:
            return lambda x: 0
        return lambda x: x ** (n - 1)  * (1 - x ** 2)


    def normalize(coefficients, f):
        p = coefficients[1] / coefficients[0]
        q = coefficients[2] / coefficients[0]
        f /= coefficients[0]
        return p, q, f
    def get_basis_system(num_of_basis_functions):
        return [Utils.get_basis_function(n) for n in range(num_of_basis_functions)]


    def get_list_of_unknown_coefficients(num):
        return [sp.Symbol('a' + str(n)) for n in range(1, num + 1)]
    def get_residual(approximate_solution, p, q, f):
        x = sp.Symbol('x')
        residual = (sp.diff(approximate_solution, x, 2) +
                    p * sp.diff(approximate_solution, x) +
                    q * approximate_solution - f)
        return residual
    def get_approximate_solution(basis_system):
        x = sp.Symbol('x')
        unknown_coefficients = Utils.get_list_of_unknown_coefficients(len(basis_system) - 1)
        approximate_solution = basis_system[0](x) + \
                               sum([unknown_coefficients[n] * basis_system[n + 1](x)\
                               for n in range(len(basis_system) - 1)])
        return approximate_solution

class Task:
    def collocation(coefficients, f, points, num_of_basis_functions):
        p, q, f = Utils.normalize(coefficients, f)
        x = sp.Symbol('x')
        basis_system = Utils.get_basis_system(num_of_basis_functions)
        unknown_coefficients = Utils.get_list_of_unknown_coefficients(len(points))
        approximate_solution = Utils.get_approximate_solution(basis_system)
        residual = Utils.get_residual(approximate_solution, p, q, f)
        system_of_equations = [residual.subs({x: point}) for point in points]
        found_coefficients = sp.solve(system_of_equations, unknown_coefficients)
        answer = approximate_solution.subs(found_coefficients)
        return answer


    def galerkin(coefficients, f, num_of_basis_functions, a, b):
        p, q, f = Utils.normalize(coefficients, f)
        x = sp.Symbol('x')
        basis_system = Utils.get_basis_system(num_of_basis_functions)
        unknown_coefficients = Utils.get_list_of_unknown_coefficients(num_of_basis_functions - 1)
        approximate_solution = Utils.get_approximate_solution(basis_system)
        residual = Utils.get_residual(approximate_solution, p, q, f)
        system_of_equations = [sp.integrate(residual * fi(x), (x, a, b)) for fi in basis_system[1:]]
        found_coefficients = sp.solve(system_of_equations, unknown_coefficients)
        answer = approximate_solution.subs(found_coefficients)
        return answer


    def integral(coefficients, f, num_of_basis_functions, a, b):
        p, q, f = Utils.normalize(coefficients, f)
        x = sp.Symbol('x')
        basis_system = Utils.get_basis_system(num_of_basis_functions)
        unknown_coefficients = Utils.get_list_of_unknown_coefficients(num_of_basis_functions - 1)
        approximate_solution = Utils.get_approximate_solution(basis_system)
        residual = Utils.get_residual(approximate_solution, p, q, f)
        I = sp.integrate(residual ** 2, (x, a, b))
        system_of_equations = [(sp.diff(I, c)) for c in unknown_coefficients]
        found_coefficients = sp.solve(system_of_equations, unknown_coefficients)
        answer = approximate_solution.subs(found_coefficients)
        return answer


    def discrete(coefficients, f, num_of_basis_functions, points):
        p, q, f = Utils.normalize(coefficients, f)
        x = sp.Symbol('x')
        basis_system = Utils.get_basis_system(num_of_basis_functions)
        unknown_coefficients = Utils.get_list_of_unknown_coefficients(num_of_basis_functions - 1)
        approximate_solution = Utils.get_approximate_solution(basis_system)
        residual = Utils.get_residual(approximate_solution, p, q, f)
        S = sum([(residual ** 2).subs({x:point}) for point in points])
        system_of_equations = [sp.diff(S, c) for c in unknown_coefficients]
        found_coefficients = sp.solve(system_of_equations, unknown_coefficients)
        answer = approximate_solution.subs(found_coefficients)
        return answer

if __name__=='__main__':
    # x = sp.Symbol('x')
    # first_equation_coefficients = [1, 0, 1 + x ** 2]
    # a = -1
    # b = 1
    # dx = 0.01
    #
    # points = list(np.linspace(-1, 1, num=10))
    # num_of_basis_functions = len(points) + 1
    #
    # answer_1_with_collocation_method = Task.collocation(
    #                                    first_equation_coefficients, -1,
    #                                    points,
    #                                    num_of_basis_functions)
    # display(answer_1_with_collocation_method)
    # Utils.plot([sp.lambdify(x, answer_1_with_collocation_method)], a, b + dx, dx)

    a = -1
    b = 1
    dx = 0.01

    x = sp.Symbol('x')
    k = 10
    coefficients = [math.sin(k), 0, 1 + math.cos(k) * x ** 2]

    points = list(np.linspace(-1, 1, num=10))
    num_of_basis_functions = len(points) + 1
    print('Collocation method')

    res_1 = Task.collocation(coefficients, -1, points, num_of_basis_functions)
    display(res_1)
    Utils.plot([sp.lambdify(x, res_1)], a, b + dx, dx)
    print('Integral method')
    res_2 = Task.integral(coefficients, -1, num_of_basis_functions, -1, 1)
    display(res_2)
    Utils.plot([sp.lambdify(x, res_2)], a, b + dx, dx)

    print('Discrete method')

    points = list(np.linspace(-1, 1, num=30))
    res_3 = Task.discrete(coefficients, -1, num_of_basis_functions, points)
    display(res_3)
    Utils.plot([sp.lambdify(x, res_3)], a, b + dx, dx)

    print('Galerkin method')

    res_4 = Task.galerkin(coefficients, -1, num_of_basis_functions, -1, 1)
    display(res_4)
    Utils.plot([sp.lambdify(x, res_4)], a, b + dx, dx)
