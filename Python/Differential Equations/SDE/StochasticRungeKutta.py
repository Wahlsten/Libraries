import numpy as np

def StochasticRungeKutta(X, a, b, b_x, dt, tp, Zp):

    Y = X + a(X, tp) * dt + b(X, tp) * np.sqrt(dt)
    Xn = X + a(X, tp) * dt + b(X, tp) * Zp + 1/2 * (b(Y, tp) - b(X, tp)) * (pow(Zp, 2) - dt) * 1./np.sqrt(dt)

    return Xn
