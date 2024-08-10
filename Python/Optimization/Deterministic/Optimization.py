from GoldenSearch           import GoldenSearch
from FibonacciOptimization  import FibonacciOptimization
from QuadraticInterpolation import QuadraticInterpolation
from NelderMead             import NelderMead
from SteepestDescent        import SteepestDescent
import numpy as np

def Optimize(interval, tol, N, method, func, func_p):
    minv = interval[0]
    maxv = interval[1]
    
    if method == 'GoldenSearch':
        [min, xmin] = GoldenSearch(interval, func, tol)
    elif method == 'FibonacciOptimization':
        [min, xmin] = FibonacciOptimization(interval, func, tol)
    elif method == 'QuadraticInterpolation':
        [min, xmin] = QuadraticInterpolation((minv + maxv)/2., func, tol)
    elif method == 'NelderMead':
        [min, xmin] = NelderMead(func, tol, N)
    elif method == 'SteepestDescent':
        [min, xmin] = SteepestDescent([minv, maxv], func, func_p, tol, N)

    return min, xmin

if __name__ == '__main__':
    #def f (x, y): return (x-2)*(x-2) + (y + 1)*(y + 1)
    def f (x): return (x-2)**4 - 1
    def fp(x): return 4*(x-2)**3

    [min, xmin] = Optimize([0, 2.1], 1e-4, 50, 'Fibonacci', f, fp)

    print('minimum =', min, 'at x =', xmin)
