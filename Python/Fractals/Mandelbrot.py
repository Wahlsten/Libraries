import numpy as np
import matplotlib.pyplot as plt
def mandelbrot(c, maxiter):
    z = c
    for n in range(maxiter):
        if abs(z) > 2:
            return maxiter - n
        z = z*z + c
    return 0

def mandelbrot_set(xmin,xmax,ymin,ymax,width,height,maxiter):
    r1 = np.linspace(xmin, xmax, width)
    r2 = np.linspace(ymin, ymax, height)
    n3 = np.empty((width,height))
    for i in range(width):
        for j in range(height):
            n3[i,j] = mandelbrot(r1[i] + 1j*r2[j], maxiter)
    return (r1,r2,n3)

xmin = -2
xmax = .6
ymin = -1.15
ymax = 1.15
x, y, z = mandelbrot_set(xmin, xmax, ymin, ymax, 1000, 1000, 100)

X, Y = np.meshgrid(x, y)
Z = np.transpose(z)

#plt.contourf(X, Y, Z, 200)
plt.imshow(-Z, extent=[xmin, xmax, ymin, ymax], cmap='RdGy')
plt.show()
