import numpy as np
import random
import matplotlib.pyplot as plt

def SierpinskiTriangle(x, y, r):
    
    if r < 1/3.:
        x1 = 0.0
        y1 = 0.0
        length = np.sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y)) / 2.0
        x_n = x1 - (x1 - x) / 2.
        y_n = y1 - (y1 - y) / 2.
    elif r < 2/3.:
        x2 = 0.5
        y2 = np.sqrt(3/4.)
        length = np.sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y)) / 2.0
        x_n = x2 - (x2 - x) / 2.
        y_n = y2 - (y2 - y) / 2.
    else:
        x3 = 1.0
        y3 = 0.0
        length = np.sqrt((x3 - x) * (x3 - x) + (y3- y) * (y3 - y)) / 2.0
        x_n = x3 - (x3 - x) / 2.
        y_n = y3 - (y3 - y) / 2.
    return x_n, y_n

def Star(x, y, r, previous):
    theta = 72/360. * 2. * np.pi

    if r < 1/5. and previous != 1:
        x1 = np.cos(0 * theta)
        y1 = np.sin(0 * theta)
        length = np.sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y)) / 2.0
        x_n = x1 - (x1 - x) / 2.
        y_n = y1 - (y1 - y) / 2.
        previous = 1
    elif r < 2/5. and previous != 2:
        x2 = np.cos(1 * theta)
        y2 = np.sin(1 * theta)
        length = np.sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y)) / 2.0
        x_n = x2 - (x2 - x) / 2.
        y_n = y2 - (y2 - y) / 2.
        previous = 2
    elif r < 3/5. and previous != 3:
        x3 = np.cos(2 * theta)
        y3 = np.sin(2 * theta)
        length = np.sqrt((x3 - x) * (x3 - x) + (y3 - y) * (y3 - y)) / 2.0
        x_n = x3 - (x3 - x) / 2.
        y_n = y3 - (y3 - y) / 2.
        previous = 3
    elif r < 4/5. and previous != 4:
        x4 = np.cos(3 * theta)
        y4 = np.sin(3 * theta)
        length = np.sqrt((x4 - x) * (x4 - x) + (y4 - y) * (y4 - y)) / 2.0
        x_n = x4 - (x4 - x) / 2.
        y_n = y4 - (y4 - y) / 2.
        previous = 4
    elif r < 5/5. and previous != 5:
        x5 = np.cos(4 * theta)
        y5 = np.sin(4 * theta)
        length = np.sqrt((x5 - x) * (x5 - x) + (y5- y) * (y5 - y)) / 2.0
        x_n = x5 - (x5 - x) / 2.
        y_n = y5 - (y5 - y) / 2.
        previous = 5
    else:
        x_n = x
        y_n = y

    return x_n, y_n, previous

def Star2(x, y, r, previous):
    theta = 72/180. * np.pi
    bias = 18/180. * np.pi
    samePrev = (previous[1] == previous[0])

    if r < 1/5. and not (((previous[0] == 2) and (previous[1] == 2)) or ((previous[0] == 5) and (previous[1] == 5))):
        x1 = np.cos(0 * theta + bias)
        y1 = np.sin(0 * theta + bias)
        length = np.sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y)) / 2.0
        x_n = x1 - (x1 - x) / 2.
        y_n = y1 - (y1 - y) / 2.
        previous_new = 1
    elif r < 2/5. and not (((previous[0] == 3) and (previous[1] == 3)) or ((previous[0] == 1) and (previous[1] == 1))):
        x2 = np.cos(1 * theta + bias)
        y2 = np.sin(1 * theta + bias)
        length = np.sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y)) / 2.0
        x_n = x2 - (x2 - x) / 2.
        y_n = y2 - (y2 - y) / 2.
        previous_new = 2
    elif r < 3/5. and not (((previous[0] == 4) and (previous[1] == 4)) or ((previous[0] == 2) and (previous[1] == 2))):
        x3 = np.cos(2 * theta + bias)
        y3 = np.sin(2 * theta + bias)
        length = np.sqrt((x3 - x) * (x3 - x) + (y3 - y) * (y3 - y)) / 2.0
        x_n = x3 - (x3 - x) / 2.
        y_n = y3 - (y3 - y) / 2.
        previous_new = 3
    elif r < 4/5. and not (((previous[0] == 5) and (previous[1] == 5)) or ((previous[0] == 3) and (previous[1] == 3))):
        x4 = np.cos(3 * theta + bias)
        y4 = np.sin(3 * theta + bias)
        length = np.sqrt((x4 - x) * (x4 - x) + (y4 - y) * (y4 - y)) / 2.0
        x_n = x4 - (x4 - x) / 2.
        y_n = y4 - (y4 - y) / 2.
        previous_new = 4
    elif r < 5/5. and not (((previous[0] == 1) and (previous[1] == 1)) or ((previous[0] == 4) and (previous[1] == 4))):
        x5 = np.cos(4 * theta + bias)
        y5 = np.sin(4 * theta + bias)
        length = np.sqrt((x5 - x) * (x5 - x) + (y5- y) * (y5 - y)) / 2.0
        x_n = x5 - (x5 - x) / 2.
        y_n = y5 - (y5 - y) / 2.
        previous_new = 5
    else:
        x_n = x
        y_n = y
        previous_new = previous[1]
        
    return x_n, y_n, previous_new

X = [0.5]
Y = [0.5]
previous = [-1, -1]

for n in range(0, 1000000):
    r = random.uniform(0, 1)
    #x, y = SierpinskiTriangle(X[n], Y[n], r)
    #x, y, previous = Star(X[n], Y[n], r, previous)
    x, y, previousTmp = Star2(X[n], Y[n], r, previous)
    previous[0] = previous[1]
    previous[1] = previousTmp
    X.append(x)
    Y.append(y)

plt.figure(figsize = [1, 1])
#plt.scatter(X, Y, color = 'b', marker = '.')
plt.plot(X, Y, linewidth = 0, color = 'k', marker = '.', markersize = 0.2)
plt.show()

