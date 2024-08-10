def LaxFriedrich(u, c):
    import numpy as np

    un = np.zeros(len(u))

    un[1:-1] = (u[:-2] +u[2:])/2.0 -  c*(u[2:] - u[:-2])/2.0

    return un[1:-1]
