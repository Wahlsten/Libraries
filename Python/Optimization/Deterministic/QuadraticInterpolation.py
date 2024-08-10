import numpy as np
def QuadraticInterpolation(x0, f, tol):

    maxj = 20
    maxk = 30
    big = 1e6
    err = 1
    k = 1
    P = []
    P.append(x0)
    cond = 0
    h = 1

    if (np.abs(x0) > 1e4):

        h = float(np.abs(x0))/1e4

    while (k < maxk) & (err > tol) & (cond != 5):

        f1 = (f(x0 + 1e-5) - f(x0 - 1e-5)) / float(2e-5)

        if (f1 > 0):
            h = -np.abs(h)

        x1 = x0 + h
        x2 = x0 + 2*h
        xmin = x0 + 0
        y0 = f(x0)
        y1 = f(x1)
        y2 = f(x2)
        ymin = y0 + 0
        cond = 0
        j = 0

        while (j < maxj) & (np.abs(h) > tol) & (cond == 0):

            if (y0 <= y1):
                x2 = x1 + 0
                y2 = y1 + 0
                h = h/2.
                x1 = x0 + h
                y1 = f(x1)
            else:
                if (y2 < y1):
                    x1 = x2 + 0
                    y1 = y2 + 0
                    h = 2*h
                    x2 = x0 + 2*h
                    y2 = f(x2)
                else:
                    cond = -1

            j = j + 1

            if (np.abs(h) > big) | (np.abs(x0) > big):
                cond = 5

        if cond == 5:
            xmin = x1 + 0
            ymin = f(x1)
        else:
            d = 4 * y1 - 2 * y0 - 2 * y2

            if (d < 0):
                hmin = h * (4 * y1 - 3 * y0 - y2)/float(d)
            else:
                hmin = h/3.
                cond = 4

            xmin = x0 + hmin
            ymin = f(xmin)
            h = np.abs(h)
            h0 = np.abs(hmin)
            h1 = np.abs(hmin - h)
            h2 = np.abs(hmin - 2*h)

            if (h0 < h):
                h = h0 + 0
            if (h1 < h):
                h = h1 + 0
            if (h2 < h):
                h = h2 + 0
            if (h == 0):
                h = hmin + 0
            if (h < tol):
                cond = 1
            if (np.abs(h) > big) | (np.abs(xmin) > big):
                cond = 5


            e0 = np.abs(x0 - ymin)
            e1 = np.abs(x1 - ymin)
            e2 = np.abs(x2 - ymin)

            if (e0 != 0) & (e0 < err):
                err = e0 + 0
            if (e1 != 0) & (e1 < err):
                err = e1 + 0
            if (e2 != 0) & (e2 < err):
                err = e2 + 0
            if (e0 != 0) & (e1 == err) & (e2 == 0):
                error = 0
            if err < tol:
                cond = 2

            x0 = xmin + 0
            k = k + 1
            P.append(x0)

        if (cond == 2) & (h < tol):
            cond = 3

    x = x0 + 0
    dx = h + 0
    dy = err + 0

    xmin = x
    minimum = f(x)

    return minimum, xmin
