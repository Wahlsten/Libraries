def DuFortFrankel(u, up, c, m):
    import numpy as np

    un = np.zeros(len(u))

    if m is 0:
        un[1:-1] = 1./(1 + 2 * c) * ((1 - 2*c) * up[1:-1] + 2*c*u[:-2] + 2 * c * u[2:])
    if m is 1:
        un[1:-1] = (1 - 2*c) * u[1:-1] + c*u[:-2] + c * u[2:]

    return un[1:-1]
