def Brent(interval, N, f, tol):
    ''' Brent '''
    import numpy as np

    a = interval[1]
    b = interval[0]
    X = []
    c = a
    d = 0
    m = True

    for n in range(1, N):

        if (f(c) != f(a)) and (f(c) != f(b)):
            s = a * f(b)*f(c)/((f(a) - f(b))*(f(a) - f(c))) + b*f(a)*f(c)/((f(b) - f(a))*(f(b) - f(c))) + c*f(a)*f(b)/((f(c) - f(a))*(f(c) - f(b)))
        else:
            s = (b * f(a) - a * f(b))/(f(a) - f(b))

        if (((s < (3.*a + b)/4.) and (s > b)) \
        or (m and (np.abs(s-b) >= np.abs(b-c)/2.)) \
        or ((not m) and (np.abs(s-b) >= np.abs(c-d)/2.)) \
        or (m and (np.abs(b-c) < tol)) \
        or ((not m) and (np.abs(c-d) < tol))):
            s = (a + b)/2.
            m = True
        else:
            m = False

        d = c
        c = b
        if f(a) * f(s) < 0:
            b = s
        else:
            a = s

        if np.abs(f(a)) < np.abs(f(b)):
            tmp = a
            a = b
            b = tmp

        X.append(s)

    return X
