import numpy as np

def GoldenSearch(interval, f, tol):

    a   = interval[0]
    b   = interval[1]
    r1  = (np.sqrt(5) - 1)/2.
    r2  = r1*r1
    h   = b-a
    ya  = f(a)
    yb  = f(b)
    c   = a + r2*h
    d   = a + r1*h
    yc  = f(c)
    yd  = f(d)

    A = []
    B = []
    C = []
    D = []
    A.append(a)
    B.append(b)
    C.append(c)
    D.append(d)

    k = 0

    while (np.abs(yb - ya) > tol) | (h > tol):
        k = k + 1

        if (yc < yd):
            b = d
            yb = yd
            d = c
            yd = yc
            h = b - a
            c = a + r2*h
            yc = f(c)
        else:
            a = c
            ya = yc
            c = d
            yc = yd
            h = b - a
            d = a + r1*h
            yd = f(d)

        A.append(a)
        B.append(b)
        C.append(c)
        D.append(d)

        dp = np.abs(b - a)
        dy = np.abs(yb - ya)
        p = a
        yp = ya

        if (yb < ya):
            p = b
            yp = yb

        G = [np.transpose(A), np.transpose(B), np.transpose(C), np.transpose(A)]
        S = [p, yp]
        E = [dp, dy]


    minimum = S[1]
    xmin    = S[0]

    return minimum, xmin
