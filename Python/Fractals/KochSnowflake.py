import numpy as np
from math import *
import matplotlib.pyplot as plt
def Koch(px, py, order):

    for n in range(len(px) - 1):

        lengthX = px[n+1] - px[n]
        lengthY = py[n+1] - py[n]
        lengthVec = sqrt(lengthX**2 + lengthY**2)

        order.append(order[n] + 1)
        order.append(order[n] + 2)
        order.append(order[n] + 3)

        px.append(px[n] + lengthX/3.0)
        py.append(py[n] + lengthY/3.0)

        new_dx = -(py[n+1] - py[n])
        new_dy = px[n+1] - px[n]

        px.append(px[n] + lengthX/2.0 + new_dx/3.0)
        py.append(py[n] + lengthY/2.0 + new_dy/3.0)

        px.append(px[n+1] - lengthX/3.0)
        py.append(py[n+1] - lengthY/3.0)

    return order, px, py

def Koch_set(px, py, levels):

    order = [k for k in range(len(px))]
    for k in range(levels):
        order = [k*4 for k in order]
        order, px, py = Koch(list(px), list(py), list(order))
        order, px, py = zip(*sorted(zip(order, px, py)))

    return px, py

px = [0, 1/2., 1, 0]
py = [0, sqrt(3/4.), 0, 0]
pxn, pyn = Koch_set(px, py, 7)

print(pxn, pyn)
plt.plot(pxn, pyn)
plt.axis('equal')
plt.show()
