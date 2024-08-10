import numpy as np
def EulerMaruyama(X, a, b, b_x, dt, tp, Zp):

    Xn = X + a(X, tp) * dt + b(X, tp) * Zp

    return Xn
