def MacCormack(u, c):
    import numpy as np

    un = np.zeros(len(u))
    up = u.copy()
    up[:-1] = u[:-1] - c*(u[1:]-u[:-1])
    un[1:] = .5*(u[1:]+up[1:] -  c*(up[1:]-up[:-1]))

    return un[1:-1]
