import numpy as np
import random
import matplotlib.pyplot as plt

def BarnsleyFern(x, y, r):
    
    if r < 1.0:
        x_n = 0
        y_n = 0.16 * y
    elif r < 86.0:
        x_n = 0.85 * x + 0.04 * y
        y_n = -0.04 * x + 0.85 * y + 1.6
    elif r < 93.0:
        x_n = 0.2 * x - 0.26 * y
        y_n = 0.23 * x + 0.22 * y + 1.6
    else:
        x_n = -0.15 * x + 0.28 * y
        y_n = 0.26 * x + 0.24 * y + 0.44

    return x_n, y_n

X = [0]
Y = [0]
for n in range(0, 100000):
    r = random.uniform(0, 100)
    x, y = BarnsleyFern(X[n], Y[n], r)
    X.append(x)
    Y.append(y)

plt.figure(figsize = [10, 10])
plt.scatter(X, Y, color = 'g', marker = '.')
plt.show()

