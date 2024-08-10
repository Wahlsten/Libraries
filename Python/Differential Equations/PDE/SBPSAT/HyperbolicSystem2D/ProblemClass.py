class ProblemClass:
    def __init__(self, prob, grid):
        import numpy as np

        if prob is 'Maxwell':
           a = np.matrix('0 0 1.; 0 0 0; 1. 0 0')
           b = np.matrix('0 -1. 0; -1. 0 0; 0 0 0')
        else:
           a = np.matrix('1. 0; 0 1.')
           b = np.matrix('1. 0; 0 1.')

        self.A = a
        self.B = b
        self.dim = a.shape[1]

    def initialData(self, data):

        def u_a(x, y, t):   return       sin(2*pi*(x - t)) + sin(2*pi*(y - t))
        def uT_a(x, y, t):  return -2*pi*cos(2*pi*(x - t)) - 2*pi*cos(2*pi*(x - t))
        def uX_a(x, y, t):  return  2*pi*cos(2*pi*(x - t))
        def uY_a(x, y, t):  return  2*pi*cos(2*pi*(y - t))

    def createData(self, grid_obj):

        def gE(x, y, t): return self.u_a(x, y, t)
        def gW(x, y, t): return self.u_a(x, y, t)
        def gN(x, y, t): return self.u_a(x, y, t)
        def gS(x, y, t): return self.u_a(x, y, t)
        def f(x, y, t):  return self.u_a(x, y, t)
