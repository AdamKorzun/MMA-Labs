from sympy import *
import numpy as np


class Task4:
    def MSI(X, Y, tol):

        x0, y0 = 1.1, 0.4
        x_guess, y_guess = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
        counter = 1
        while (abs(x0 - x_guess) > tol or abs(y0 - y_guess) > tol):
            x0, y0 = x_guess, y_guess
            x_guess, y_guess = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [x_guess, y_guess]

    def NewtonM(f1, f2,x, y, tol):
        matrix = Matrix([f1, f2])
        J = Matrix([[diff(f1, x), diff(f1, y)], [diff(f2, x), diff(f2, y)]]).inv()
        est = Matrix([1.1, 0.4])
        guess = (est - J.inv()*matrix).subs(x, est[0]).subs(y, est[1])
        counter = 1
        while (abs(est[0] - guess[0]) > tol):
            est = guess
            guess = (est - J*matrix).subs(x, est[0]).subs(y, est[1])
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [guess[0], guess[1]]
def test_1(tol):
    x = Symbol('x')
    y = Symbol('y')
    func1 = x**2 - 2*y**2 - x * y + 2*x - y + 1
    func2 = 2*x**2 - y**2 + x * y + 3*y - 5
    sol = nsolve([func1, func2],[x,y], [-1,-1])
    print(func1.subs([(x,sol[0]),(y, sol[1])]))
    X = solve(func1, x)[1]
    Y = solve(func2, y)[1]
    '''
    p1 = plot(solve(func1, y)[0], show = False)
    p2 = plot(solve(func2, y)[0], show = False)
    p3 = plot(solve(func1, y)[1], show = False)
    p4 = plot(solve(func2, y)[1], show = False)
    p1.append(p2[0])
    p1.append(p3[0])
    p1.append(p4[0])
    p1.show()
    '''
    print(f'{X}\n{Y}')
    print(Task4.NewtonM(func1, func2, x, y, tol))
    print(Task4.MSI(X, Y, tol))
def test_2(tol):
        x = Symbol('x')
        y = Symbol('y')
        func1 = 3*x**2 + x*y - 2 * y**2
        func2 = 2 * x**2 - 3*y + y**2 + 1
        X = solve(func1, x)[1]
        Y = solve(func2, y)[1]
        print(f'{X}\n{Y}')
        sol = Task4.NewtonM(func1, func2, x, y, tol)
        print(func1.subs([(x,sol[0]),(y, sol[1])]))
        print(Task4.MSI(X, Y, tol))
def test_3(tol):
        x = Symbol('x')
        y = Symbol('y')
        
        sol = nsolve([func1, func2],[x,y], [-1,-1])
        print(func1.subs([(x,sol[0]),(y, sol[1])]))
        X = solve(func1, x)[0]
        Y = solve(func2, y)[0]
        print(f'{X}\n{Y}')
        sol = Task4.NewtonM(func1, func2, x, y, tol)
        print(func1.subs([(x,sol[0]),(y, sol[1])]))
        print(Task4.MSI(X, Y, tol))
if __name__ == "__main__":
    tol = 0.0001
    a = 0.6
    m = 0.4
    x = Symbol('x')
    y = Symbol('y')
    func1 = tan(x*y + m) - x
    func2 = a*x**2 + 2*y**2 - 1



    X = sqrt((1-2*y**2)/a)
    Y = (atan(x) - m) / x
    #print(test_1(tol))
    #print(test_2(tol))
    print(test_3(tol))
    print(Task4.MSI(X, Y, tol))
    print(Task4.NewtonM(func1, func2, x, y, tol))
    p1 = plot(sqrt((1-a*x**2)/2), show = False)
    p2 = plot((atan(x) - m) / x, show = False)
    p1.append(p2[0])
    p1.show()
