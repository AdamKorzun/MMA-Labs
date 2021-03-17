from sympy import *
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
        x1, y1 = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
        counter = 1
        while (abs(x0 - x1) > tol or abs(y0 - y1) > tol):
            x0, y0 = x1, y1
            x1, y1 = X.subs(x, x0).subs(y, y0), Y.subs(x, x0).subs(y, y0)
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [x1, y1]
    def NewtonM(self, tol):
        a = self.a
        m = self.m
        x = self.x
        y = self.y
        X = self.X
        Y = self.Y
        f1 = tan(x*y + m) - x
        f2 = a*x**2 + 2*y**2 - 1
        F = Matrix([f1, f2])
        J = Matrix([[diff(f1, x), diff(f1, y)], [diff(f2, x), diff(f2, y)]])
        v0 = Matrix([1.1, 0.4])
        J = J.inv()
        v1 = (v0 - J.inv()*F).subs(x, v0[0]).subs(y, v0[1])
        counter = 1
        while (abs(v0[0] - v1[0]) > tol or abs(v0[1] - v1[1]) > tol):
            v0 = v1
            v1 = (v0 - J*F).subs(x, v0[0]).subs(y, v0[1])
            counter += 1
        print('Number of iterations: {}'.format(counter))
        return [v1[0], v1[1]]

if __name__ == "__main__":
    tol = 0.0001
    a = 0.6
    m = 0.4
    x = Symbol('x')
    y = Symbol('y')
    task4 = Task4(a, m, x, y, sqrt((1-2*y**2)/a), (atan(x) - m) / x)

    print(task4.MSI(tol))
    print(task4.NewtonM(tol))
    p1 = plot(sqrt((1-2*y**2)/a), show = False)
    p2 = plot((atan(x) - m) / x, show = False)
    p1.append(p2[0])
    p1.show()
