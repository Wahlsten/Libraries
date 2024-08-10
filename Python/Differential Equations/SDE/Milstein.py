import numpy as np

def Milstein(X, a, b, b_x, dt, tp, Zp):

    Xn = X + a(X, tp) * dt + b(X, tp) * Zp + 1/2 * b_x(X, tp) * (pow(Zp, 2) - dt)

    return Xn
