def NewtonRaphson(interval, N, f, fp):
    ''' Newton Raphson '''
    import numpy as np

    xp = (interval[1] + interval[0])/2
    X = []

    for n in range(1, N):

        xp = xp - f(xp)/fp(xp)
        X.append(xp)

    return X
