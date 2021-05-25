import numpy as np
class Types:
    simpson = 'simpson'
    right_rect = 'right_rect'
    left_rect = 'left_rect'
    trapezium = 'trapezium'
    central_rect = 'central_rect'

class Task8:


    def it(lb, rb, tol, power):
        return (rb - lb) / np.power(tol, power)

    def calc_integral(f, left_border, right_border, tol, option):
        if option == Types.simpson:
            n = np.ceil(Task8.it(left_border, right_border, tol, 0.25))
            h = (right_border - left_border) / n
            integral = Task8.integration_types[option](f, n, left_border, h)
            previous_integral = integral + tol * 20
            counter = 0
            while abs(integral - previous_integral)  > tol:
                counter += 1
                previous_integral = integral
                h /= 2
                n *= 2
                integral = Task8.integration_types[option](f, n, left_border, h)
            return integral, counter
        else:
            n = np.ceil((right_border - left_border) / np.power(tol, 1 / 4))
            h = (right_border - left_border) / n
            integral = Task8.integration_types[option](f, n, left_border, h)
            previous_integral = integral + tol
            counter = 0
            while True:
                counter += 1
                if abs(integral - previous_integral) < 0.0000001:
                    return integral, counter
                previous_integral = integral
                h /= 2
                n *= 2
                integral = Task8.integration_types[option](f, n, left_border, h)
            return integral, counter

    integration_types = dict(
        left_rect=lambda f, n, left_border, h:
        sum([f(left_border + i * h) for i in range(int(n))]) * h,
        right_rect=lambda f, n, left_border, h:
        sum([f(left_border + (i + 1) * h) for i in range(int(n))]) * h,
        central_rect=lambda f, n, left_border, h:
        sum([f(left_border + (i + 1 / 2) * h) for i in range(int(n))]) * h,
        trapezium=lambda f, n, left_border, h:
        (sum([f(left_border + i * h) for i in range(1, int(n))]) + 0.5 *
         (f(left_border + h * n) + f(left_border))) * h,
        simpson=lambda f, n, left_border, h:
        (h / 3) * (f(left_border) + f(left_border + float(n) * h)
                   + 2 * sum([f(2 * float(i) * h) for i in range(1, int(n) // 2)])
                   + 4 * sum([f(2 * i + 1) for i in range(int(n) // 2)])))

    def first_diff(f, x, tol):
        diff = 10 * tol
        h = 0.1
        derivative = (f(x + h) - f(x - h)) / (2 * h)
        counter = 0
        while abs(diff) > tol:
            temp = derivative
            derivative = (f(x + h) - f(x - h)) / (2 * h)
            if counter != 0:
                diff = derivative - temp
            else:
                diff = 10 * tol
            h /= 4
            counter += 1
        return derivative, counter

    def second_diff(f, x, tol):
        counter = 0
        h = 0.1
        diff = 10 * tol
        derivative = (f(x + h) - f(x - h)) / (2 * h)

        while abs(diff) > tol:
            temp = derivative
            derivative = (f(x + h) + f(x - h) - 2 * f(x)) / (4 * h ** 2)
            if counter != 0:
                diff = derivative - temp
            else:
                diff = 10 * tol
            h /= 4
            counter += 1
        return derivative, counter


def t():
    f = lambda x: np.sin(x) / x
    print(Task8.calc_integral(f, 1, 2, 0.000001, Types.left_rect))
    print(Task8.calc_integral(f, 1, 2, 0.000001, Types.right_rect))
    print(Task8.calc_integral(f, 1, 2, 0.000001, Types.central_rect))
    print(Task8.calc_integral(f, 1, 2, 0.000001, Types.trapezium))
    print(Task8.calc_integral(f, 1, 2, 0.000001, Types.simpson)).
    
if __name__ == "__main__":
    f = lambda x: np.sin(x) / x
    t()
    print(Task8.first_diff(f, 1.5, 0.01))
    print(Task8.second_diff(f, 1.5, 0.01))
