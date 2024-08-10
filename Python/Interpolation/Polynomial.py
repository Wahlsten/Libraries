import numpy as np
def Polynomial(x0, y0, xn):
    print('Polynomial')

    N0 = len(x0)
    NN = len(xn)
    yn = np.zeros(NN)
    h = (x0[1:] - x0[0:-1])
    k = 0

    A = np.zeros((N0, N0))

    for j in range(0, N0):

        A[:,j] = np.transpose(x0**j)

    a = np.linalg.solve(A,y0)

    for k in range(0, N0):

        yn = yn + a[k]*xn**k

    return yn
