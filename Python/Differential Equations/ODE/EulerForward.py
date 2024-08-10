import numpy as np
def EulerForward(interval, f0, ft, N):
    h = (interval[1] - interval[0])/float(N)
    y = np.zeros((N + 1, 1))
    y = y[:,0]
    t = np.linspace(interval[0], interval[1], N + 1, endpoint = True)
    y[0] = f0

    for j in range(0, N):
        y[j+1] = y[j] + h*ft(t[j], y[j])

    return t, y
