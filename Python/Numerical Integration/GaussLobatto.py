def GaussLobatto(interval, f, N):
    ''' Gauss-Lobatto '''
    import numpy as np

    if N == 3:
        p = np.array([0, 1, -1])
        w = np.array([4/3., 1/3., 1/3.])
    elif N == 4:
        p = np.array([np.sqrt(1/5.), -np.sqrt(1/5.), 1, -1])
        w = np.array([5/6., 5/6., 1/6., 1/6.])
    elif N == 5:
        p = np.array([0, np.sqrt(3./7.), -np.sqrt(3./7.), 1, -1])
        w = np.array([32/45., 49/90., 49/90., 1/10., 1/10.])
    elif N == 6:
        p = np.array([np.sqrt(1/3. - 2*np.sqrt(7)/21.), -np.sqrt(1/3. - 2*np.sqrt(7)/21.), np.sqrt(1/3. + 2*np.sqrt(7)/21.), -np.sqrt(1/3. + 2*np.sqrt(7)/21.), 1, -1])
        w = np.array([(14 + np.sqrt(7))/30., (14 + np.sqrt(7))/30., (14 - np.sqrt(7))/30., (14 - np.sqrt(7))/30., 1/15., 1/15.])
    elif N == 7:
        p = np.array([0, np.sqrt(5/11. - 2/11.*np.sqrt(5/3.)), -np.sqrt(5/11. - 2/11.*np.sqrt(5/3.)), np.sqrt(5/11. + 2/11.*np.sqrt(5/3.)), -np.sqrt(5/11. + 2/11.*np.sqrt(5/3.)), 1, -1])
        w = np.array([256/525., (124 + 7*np.sqrt(15))/350., (124 + 7*np.sqrt(15))/350., (124 - 7*np.sqrt(15))/350., (124 - 7*np.sqrt(15))/350., 1/21., 1/21.])
    else:
        print('number of points undefined')

    a = (interval[1] - interval[0])/2.
    b = (interval[1] + interval[0])/2.

    I = a * np.dot(f(a*p + b), w)

    return I
