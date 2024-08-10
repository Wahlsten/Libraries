def CentralTimeCentralSpace(u, up, c):
    import numpy as np

    un = np.zeros(len(u))

    un[1:-1] = (2 - 2*c) * u[1:-1] - up[1:-1] + c*u[:-2] + c * u[2:]

    return un[1:-1]
