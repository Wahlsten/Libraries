def Downwind(u, c, d):
    import numpy as np

    un = np.zeros(len(u))

    un[1:-1] = (1 - 2 * d + c) * u[1:-1] + d * u[:-2] + (-c + d) * u[2:]

    return un[1:-1]
