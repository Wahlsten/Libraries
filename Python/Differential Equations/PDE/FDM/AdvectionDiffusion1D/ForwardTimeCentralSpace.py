def ForwardTimeCentralSpace(u, c, d):
    import numpy as np

    un = np.zeros(len(u))

    un[1:-1] = (1 - 2 * d) * u[1:-1] + (c/2. + d) * u[:-2] + (-c/2. + d) * u[2:]

    return un[1:-1]
