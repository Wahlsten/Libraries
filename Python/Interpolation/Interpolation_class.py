from Constant   import Constant
from Linear     import Linear
from Polynomial import Polynomial
from Splines    import Splines        
import matplotlib.pyplot as plt
import numpy as np
from math import *

def Interpolate(x0, y0, xn, method):

    if method == 'Constant':
        yn = Constant(x0, y0, xn)
    elif method == 'Linear':
        yn = Linear(x0, y0, xn)
    elif method == 'Polynomial':
        yn = Polynomial(x0, y0, xn)
    elif method == 'Splines':
        yn = Splines(x0, y0, xn)
    return yn

def f(x): return np.sin(pi*x)

x0 = np.linspace(-4,4,20, endpoint = True)
y0 = f(x0)
xn = np.linspace(-4,4,100, endpoint = True)

yn = Interpolate(x0, y0, xn, 'Polynomial')

plt.plot(xn, yn)
plt.show()

plt.plot(xn, f(xn))
plt.show()

plt.plot(xn, f(xn), xn, yn)
plt.show()
