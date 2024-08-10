import numpy as np

def Bisection(interval, N, f):
    ''' Bisection '''

    upper = interval[1]
    lower = interval[0]
    X = []

    for n in range(1, N):

        middle = (upper + lower)/2

        if f(middle) * f(upper) > 0:
            upper = middle
        else:
            lower = middle

        X.append(middle)

    return X
