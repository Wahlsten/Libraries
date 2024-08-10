import numpy as np

def SteepestDescent(interval, f, fp, tol, N):

    def g(x): return -(1/np.linalg.norm(fp(x))) * fp(x)

    maxj = 10
    big = 1e8
    h = 1
    n = 1
    P = np.zeros((N, n+1))
    x0 = (interval[1] - interval[0])/2.
    leng = np.linalg.norm(x0)
    y0 = f(x0)
    delta = tol

    if (leng > 1e4):
        h = leng/float(1e4)

    err = 1
    cond = 0
    P[0, :] = [x0, y0]
    cnt = 1
    while ((fp(x0) != 0) & (cnt < N) & (cond != 5)) & ((h > delta) | (err > tol)):

        S = g(x0)
        x1 = x0 + h*S
        x2 = x0 + 2*h*S
        y1 = f(x1)
        y2 = f(x2)
        cond = 0
        j = 0

        while (j < maxj) & (cond == 0):
            leng = np.linalg.norm(x0)
            if (y0 < y1):
                x2 = x1+0
                y2 = y1+0
                h = h/2.
                x1 = x0 + h*S
                y1 = f(x1)
            else:
                if (y2 < y1):
                    x1 = x2+0
                    y1 = y2+0
                    h = 2.*h
                    x2 = x0 + 2*h*S
                    y2 = f(x2)
                else:
                    cond = -1

            j = j+1
            if (h < delta):
                cond = 1

            if (np.abs(h) > big) | (leng > big):
                cond = 5

        if (cond == 5):
            xmin = x1+0
            ymin = y1+0

        else:
            d = 4*y1 - 2*y0 - 2*y2
            if (d < 0):
                hmin = h * (4 * y1 - 3 * y0 - y2)/d
            else:
                cond = 4
                hmin = h/3

            xmin = x0 + hmin*S
            ymin = f(xmin)
            h0 = np.abs(hmin)
            h1 = np.abs(hmin - h)
            h2 = np.abs(hmin - 2 * h)

            if (h0 < h):
                h = h0+0
            if (h1 < h):
                h = h1+0
            if (h2 < h):
                h = h2+0
            if (h == 0):
                h = hmin+0
            if (h < delta):
                cond = 1

            e0 = np.abs(y0 - ymin)
            e1 = np.abs(y1 - ymin)
            e2 = np.abs(y2 - ymin)

            if (e0 != 0) & (e0 < err):
                err = e0+0
            if (e1 != 0) & (e1 < err):
                err = e1+0
            if (e2 != 0) & (e2 < err):
                err = e2+0
            if (e0 == 0) & (e1 == 0) & (e2 == 0):
                err = 0
            if (err < tol):
                cond = 2
            if (cond == 2) & (h < delta):
                cond = 3

        P[cnt, :] = [xmin, ymin]
        cnt = cnt + 1
        x0 = xmin
        y0 = ymin


    minimum  = ymin
    return minimum, xmin
