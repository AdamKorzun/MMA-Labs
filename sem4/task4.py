from sympy import *
import numpy as np


class Task4:
    def __init__(self, a, m,x, y, X, Y):
        self.a = a
        self.m = m
        self.x = x
        self.y = y
        self.X = X
        self.Y = Y
    def MSI(self, tol):
        a = self.a
        m = self.m
        x = self.x
        y = self.y
        X = self.X
        Y = self.Y
        x0, y0 = 1.1, 0.4
        x_guess, y_guess = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
        counter = 1
        while (abs(x0 - x_guess) > tol or abs(y0 - y_guess) > tol):
            x0, y0 = x_guess, y_guess
            x_guess, y_guess = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [x_guess, y_guess]
    def NewtonM(self,f1, f2, tol):
        a = self.a
        m = self.m
        x = self.x
        y = self.y
        F = Matrix([f1, f2])
        J = Matrix([[diff(f1, x), diff(f1, y)], [diff(f2, x), diff(f2, y)]])
        est = Matrix([1.1, 0.4])
        J = J.inv()
        guess = (est - J.inv()*F).subs(x, est[0]).subs(y, est[1])
        counter = 1
        #or abs(est[1] - guess[1]) > tol
        while (abs(est[0] - guess[0]) > tol):
            est = guess
            guess = (est - J*F).subs(x, est[0]).subs(y, est[1])
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [guess[0], guess[1]]

if __name__ == "__main__":
    tol = 0.0001
    a = 0.6
    m = 0.4
    x = Symbol('x')
    y = Symbol('y')
    task4 = Task4(a, m, x, y, sqrt((1-2*y**2)/a), (atan(x) - m) / x)

    print(task4.MSI(tol))
    print(task4.NewtonM(tan(x*y + m) - x, a*x**2 + 2*y**2 - 1, tol))
    p1 = plot(sqrt((1-2*y**2)/a), show = False)
    p2 = plot((atan(x) - m) / x, show = False)
    p1.append(p2[0])
    p1.show()
