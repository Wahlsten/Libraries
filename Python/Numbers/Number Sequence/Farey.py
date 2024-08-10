import numpy as np
def Farey(N):
    ''' Farey '''

    a, b, c, d = 0, 1, 1, N
    Num = []
    Num.append(float(a)/b)
    while (c <= N):

        k = int((N + b) / d)
        a, b, c, d = c, d, (k*c-a), (k*d-b)
        Num.append(float(a)/b)

    return Num
