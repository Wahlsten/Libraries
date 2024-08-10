import numpy as np
import matplotlib.pyplot as plt

def JuliaSet(z_start, c, max_iterations):
    z = z_start
    for iterations in range(max_iterations):
        if abs(z) > 2:
            return max_iterations - iterations
        z = z**2 + c
        #z = z*z + c
        #z = c * np.sin(z) #or z = c*cos(z)

    return 0

def Julia_set(xmin, xmax, ymin, ymax, width, height, maxiter):
    x_vals = np.linspace(xmin, xmax, width)
    y_vals = np.linspace(ymin, ymax, height)
    #data   = np.zeros((width, height), dtype=np.uint8)
    data   = np.zeros((width, height))
    c      = -0.4 + 0.6j
    for x_index, x in enumerate(x_vals):
        for y_index, y in enumerate(y_vals):
            z = complex(x, y)
            data[x_index, y_index] = JuliaSet(z, c, maxiter)

    return (x, y, data)

xmin = -2
xmax =  2
ymin = -2
ymax =  2
x, y, z = Julia_set(xmin, xmax, ymin, ymax, 2000, 2000, 100)

X, Y = np.meshgrid(x, y)
Z = np.transpose(z)

#plt.contourf(X, Y, Z, 200)
plt.imshow(-Z, extent=[xmin, xmax, ymin, ymax], cmap='coolwarm')
plt.show()
