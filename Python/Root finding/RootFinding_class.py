from Bisection import Bisection
from RegulaFalsi import RegulaFalsi
from Secant import Secant
from NewtonRaphson import NewtonRaphson
from Brent import Brent

def Solve(minv, maxv, N, method, f, fp, tol):

    minv = float(minv)
    maxv = float(maxv)
    
    if method == 'Bisection':
        x_n = Bisection([minv, maxv], N, f)
    elif method == 'RegulaFalsi':
        x_n = RegulaFalsi([minv, maxv], N, f)
    elif method == 'Secant':
        x_n = Secant([minv, maxv], N, f)
    elif method == 'NewtonRaphson':
        x_n = NewtonRaphson([minv, maxv], N, f, fp)
    elif method == 'Brent':
        x_n = Brent([minv, maxv], N, f, tol)
    return x_n

if __name__ == '__main__':
    def f (x): return x**3 - 27
    def fp(x): return 3*x**2

    x_n = Solve(1, 10, 10, 'Brent', f, fp, 1e-4)

    print(x_n)
