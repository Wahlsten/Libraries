def MidpointRule(interval, f, N):
    ''' Midpoint rule '''
    import numpy as np

    h = (interval[1] - interval[0])/float(N)

    grid = np.linspace(interval[0] + h/2., interval[1] - h/2., N, endpoint=True)

    I = h * (np.sum(f(grid[0:])))

    return I
