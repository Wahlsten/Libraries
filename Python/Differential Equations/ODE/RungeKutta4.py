import numpy as np
def RungeKutta4(interval, f0, ft, N):
    dim = len(interval) - 1
    h = (interval[1] - interval[0])/float(N)
    y = np.array(np.zeros((dim, N + 1)))
    t = np.linspace(interval[0], interval[1], N + 1, endpoint = True)
    y[:, 0] = f0

    for j in range(0, N):

        k1 = ft(t[j], y[:, j])
        k2 = ft(t[j] + h/2., y[:, j] + h / 2. * k1)
        k3 = ft(t[j] + h/2., y[:, j] + h / 2. * k2)
        k4 = ft(t[j] + h, y[:, j] + h * k3)
        y[:, j + 1] = y[:, j] + h / 6.*(k1 + 2 * k2 + 2 * k3 + k4)

    return t, y[0]
