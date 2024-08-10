import numpy as np

def FibonacciOptimization(interval, f, tol):

    def fib(x):

        if x < 3:
            return 1
        else:
            return fib(x - 1) + fib(x - 2)

    a = interval[0]
    b = interval[1]
    e = 0
    i = 1
    F = 1

    while F <= (b-a)/tol:
        F = fib(i)
        i = i + 1

    n = i - 1
    A = np.zeros((n - 2, 1))
    B = np.zeros((n - 2, 1))
    A[0] = a
    B[0] = b
    c = A[0] + (float(fib(n-2))/fib(n)) * (B[0] - A[0])
    d = A[0] + (float(fib(n-1))/fib(n)) * (B[0] - A[0])
    k = 0

    while k < n - 3:

        if f(c) > f(d):
            A[k + 1] = c
            B[k + 1] = B[k]
            c = d
            d = A[k+1] + (float(fib(n - k - 1))/fib(n - k)) * (B[k+1] - A[k+1])

        else:
            A[k+1] = A[k]
            B[k+1] = d
            d = c
            c = A[k+1] + (float(fib(n-k-2))/fib(n-k)) * (B[k+1] - A[k+1])

        k = k + 1

    if f(c) > f(d):
        A[n-3] = c
        B[n-3] = B[n-4]
        c = d
        d = A[n-3] + (0.5 + e) * (B[n-3] - A[n-3])
    else:
        A[n-3] = A[n-4]
        B[n-3] = d
        d = c
        c = A[n-3] + (0.5 - e) * (B[n-3] - A[n-3])

    if f(c) > f(d):
        a = c
        b = B[n-3]
    else:
        a = A[n-3]
        b = d

    minimum = f((a[0] + b[0])/2)
    xmin    = (a[0] + b[0])/2

    return minimum, xmin
