import numpy as np

def Splines(x0, y0, xn):

    N0 = len(x0)
    NN = len(xn)
    yn = np.zeros(NN)
    a  = np.zeros(N0)
    b  = np.zeros(N0)
    c  = np.zeros(N0)
    d  = np.zeros(N0)
    e  = np.zeros(N0)
    f  = np.zeros(N0)
    g  = np.zeros(N0)

    h = (x0[1:] - x0[0:-1])
    a = 3./h * (y0[1:] - y0[:-1])
    a[1:] = a[1:] - 3./h[:-1] * (y0[1:-1] - y0[:-2])

    e[0] = 1
    f[0] = 0
    g[0] = 0

    for k in range(1, N0-1):

        e[k] = 2*(x0[k+1] - x0[k-1]) - h[k-1]*f[k-1]
        f[k] = h[k]/e[k]
        g[k] = (a[k] - h[k-1]*g[k-1])/e[k]

    e[-1] = 1
    g[-1] = 0
    c[-1] = 0

    for l in range(N0-2, -1, -1):

        c[l] = g[l] - f[l] * c[l+1]
        b[l] = (y0[l+1] - y0[l])/h[l] - h[l]*(c[l+1] + 2*c[l])/3.
        d[l] = (c[l+1] - c[l])/(3.*h[l])

    S = np.array([y0, b, c, d, x0])
    S = np.transpose(S)

    m = 0

    for j in range(0, NN):

        for i in range(m, N0-1):

            if (xn[j] >= x0[i]) & (xn[j] <= x0[i+1]):
                yn[j] = y0[i] + b[i]*(xn[j] - x0[i]) + c[i]*(xn[j] - x0[i])**2 + d[i]*(xn[j] - x0[i])**3
                m = i + 0
                break

    return yn
