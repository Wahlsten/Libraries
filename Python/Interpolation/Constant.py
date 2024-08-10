import numpy as np
def Constant(x0, y0, xn):
    print('Constant')

    N0 = len(x0)
    NN = len(xn)
    yn = np.zeros(NN)
    h = np.zeros(N0+1)
    h[0] = x0[1] - x0[0]
    h[1:-1] = (x0[1:] - x0[0:-1])
    h[-1] = h[-2] + 0
    k = 0

    for j in range(0, NN):

        for i in range(k, N0):

            if (xn[j] > x0[i] - h[i]/2.) & (xn[j] < x0[i] + h[i+1]/2.):
                yn[j] = y0[i] + 0
                k = i + 0
                break

    return yn
