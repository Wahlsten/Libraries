def GaussLegendre(interval, f, N):
    ''' Gauss-Legendre '''
    import numpy as np

    if N == 1:
        p = np.array([0])
        w = np.array([2])
    elif N == 2:
        p = np.array([0.57735, -.57735])
        w = np.array([1, 1])
    elif N == 3:
        p = np.array([0, 0.774597, -0.774597])
        w = np.array([0.888889, 0.555556, 0.555556])
    elif N == 4:
        p = np.array([0.339981, -0.339981, 0.861136, -0.861136])
        w = np.array([0.652145, 0.652145, 0.347855, 0.347855])
    elif N == 5:
        p = np.array([0, 0.538469, -0.538469, 0.90618, -0.90618])
        w = np.array([0.568889, 0.478629, 0.478629, 0.236927, 0.236927])
    else:
        print('number of points undefined')

    a = (interval[1] - interval[0])/2.
    b = (interval[1] + interval[0])/2.

    I = a * np.dot(f(a*p + b), w)

    return I
