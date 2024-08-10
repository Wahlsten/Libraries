import numpy as np

def Linear(x0, y0, xn):
    print('Linear')

    N0 = len(x0)
    NN = len(xn)
    yn = np.zeros(NN)
    h = (x0[1:] - x0[0:-1])
    k = 0

    for j in range(0, NN):

        for i in range(k, N0-1):

            if (xn[j] >= x0[i]) & (xn[j] <= x0[i+1]):
                yn[j] = y0[i] + (y0[i+1]  - y0[i])*(xn[j] - x0[i])/(x0[i+1] - x0[i])
                k = i + 0
                break

    return yn
