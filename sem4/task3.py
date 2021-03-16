from sympy import sturm
from sympy.abc import x
from sympy import Poly
import numpy as np
from scipy import misc

def sturm_theorem(polyn, root_range):
    if len(root_range) != 2:
        raise Exception('invalid root range')
    sturm_seq = sturm(polyn)
    sign_left = 0
    for index, ran_num in enumerate(root_range):
        subs_seq = []
        sign_changes = 0
        for el in sturm_seq:
            subs_num = el.subs({x:ran_num})
            subs_seq.append(subs_num)
        subs_seq = np.array(subs_seq)
        for i in range(1, len(subs_seq)):
            if np.sign(subs_seq[i - 1]) != np.sign(subs_seq[i]):
                sign_changes+=1
        if index == 1:
            return sign_left - sign_changes
        else:
            sign_left = sign_changes


def separate_roots(polynom, root_range):
    number_of_roots = sturm_theorem(polynom, root_range)
    mid  = (root_range[1] - root_range[0]) / 2
    if number_of_roots == 1:
        return root_range
    if np.abs(mid) < 0.01:
        return
    roots = separate_roots(polynom, [root_range[0] + mid, root_range[1]]) # needs fix
    roots_2 = separate_roots(polynom, [root_range[0],root_range[1] - mid])
    if roots is not None:
        return roots
    elif roots_2 is not None:
        return roots_2
    return

f = lambda x: x**3 + -13.3667 * x**2 +39.8645 * x + -20.6282

def bisection(f, a, b, tol, counter = 0):
    counter +=1
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("no roots found")
    m = (a + b)/2
    if np.abs(f(m)) < tol:
        print('Number of iterations: {}'.format(counter))
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        return bisection(f, m, b, tol, counter)
    elif np.sign(f(b)) == np.sign(f(m)):
        return bisection(f, a, m, tol, counter)

def secant(f,a,b,tol):
    if f(a)*f(b) >= 0:
        raise Exception('Error')
    a_n = a
    b_n = b
    max_iterations = 1000
    for n in range(1,max_iterations):
        m_n = a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))
        f_m_n = f(m_n)
        if np.abs(f_m_n)  <  tol:
            print('Number of iterations: {}'.format(n))
            return m_n
        elif f(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        else:
            raise Exception('Error')
    return a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))

def newtons_method(f, x, tol):
    while True:
        x1 = x - f(x) / misc.derivative(f, x)
        t = abs(x1 - x)
        if t < tolerance:
            break
        x = x1
    return x

if __name__ == "__main__":
    a = -13.3667
    b = 39.8645
    c = -20.6282
    root_range = [-10, 10]
    polynom = Poly(x**3 + a * x**2 +b * x + c)
    number_of_roots = sturm_theorem(polynom, root_range)
    ranges = []
    for i in range(number_of_roots):
        sep_range = separate_roots(polynom, root_range.copy())
        print(root_range)
        print(sep_range)

        if sep_range[0] == root_range[0]:
            root_range[0] = sep_range[1]
        else:
            root_range[1] = sep_range[0]
        ranges.append(sep_range.copy())





    print(ranges)
    tol = 0.00001
    print(bisection(f, ranges[0][0], ranges[0][1], tol))
    print(secant(f, ranges[0][0], ranges[0][1], tol))
    print(newtons_method(f, -10, tol))
