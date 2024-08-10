def LaxWendroff(u, c):
    import numpy as np

    un = np.zeros(len(u))

    un[1:-1] = c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]

    return un[1:-1]
