import numpy as np

def SimpsonsRule(interval, f, N):
    ''' Simpsons rule '''

    h = (interval[1] - interval[0])/float(N)

    grid = np.linspace(interval[0], interval[1], N + 1, endpoint=True)

    I = 2.*h/6. * (f(grid[0]) + 4*np.sum(f(grid[1:-1:2])) + 2*np.sum(f(grid[2:-1:2])) + f(grid[-1]))

    return I
