def Romberg(interval, f, N, tol):
    ''' Romberg integration '''
    import numpy as np
    import math

    M = 1
    h = interval[1] - interval[0]
    err = 1
    J = -1
    R = np.zeros((4,4))
    R[0, 0] = h * (f(interval[1]) + f(interval[0])) / 2.

    while ((err > tol) & (J < 2)) | (J < 2):
        J = J + 1
        h = h/2.
        s = 0

        for p in range(1, M+1):

            x = interval[0] + h*(2*p - 1)
            s = s + f(x)

        R[J+1, 0] = R[J, 0]/2. + h*s
        M = 2*M

        for K in range(0, J+1):
            R[J+1, K+1] = R[J+1, K] + (R[J+1, K] - R[J, K])/(math.pow(4.0, K+1) - 1)

        err = np.abs(R[J, J] - R[J+1, J+1])

    I = R[J+1, J+1]

    return I
