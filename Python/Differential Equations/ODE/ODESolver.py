import numpy as np
import matplotlib.pyplot as plt
from Heun         import Heun
from RungeKutta4  import RungeKutta4
from EulerForward import EulerForward

def ODESolve(minv, maxv, N, method, ft, f0):

    if method == 'EulerForward':
        [t, y] = EulerForward([minv, maxv], f0, ft, N)
    elif method == 'Heun':        
        [t, y] = Heun([minv, maxv], f0, ft, N)
    elif method == 'RungeKutta4':        
        [t, y] = RungeKutta4([minv, maxv], f0, ft, N)

    return t, y

if __name__ == '__main__':
    def y_a(t):   return np.sin(t)
    def ft(t, y): return y - np.sin(t) + np.cos(t)

    [t, y] = ODESolve(0, 1, 100, 'RungeKutta4', ft, y_a(0))
    plt.plot(t, y)
    plt.plot(t, y_a(t))
    plt.show()
    print(abs(y[-1] - y_a(1)))
