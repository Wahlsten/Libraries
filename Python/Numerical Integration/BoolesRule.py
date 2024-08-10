def BoolesRule(interval, f, N):
    ''' Booles rule '''
    import numpy as np

    h = (interval[1] - interval[0])/float(N)

    grid = np.linspace(interval[0], interval[1], N + 1, endpoint=True)

    I = 2.*h/45. * (7*f(grid[0]) + 32*np.sum(f(grid[1:-1:4])) \
    + 12*np.sum(f(grid[2:-1:4])) + 32*np.sum(f(grid[3:-1:4])) \
    + 14*np.sum(f(grid[4:-1:4])) + 7*f(grid[-1]))

    return I
