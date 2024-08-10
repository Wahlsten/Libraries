def RegulaFalsi(interval, N, f):
    ''' Regula Falsi '''
    import numpy as np

    upper = interval[1]
    lower = interval[0]
    X = []

    for n in range(1, N):

        middle = (lower * f(upper) - upper * f(lower))/(f(upper) - f(lower))

        if f(middle) * f(upper) > 0:
            upper = middle
        else:
            lower = middle

        X.append(middle)

    return X
