class ProblemClass:
    def __init__(self, prob, grid):
        import numpy as np

        if prob is 'Maxwell':
           a = np.matrix('0 1.; 1. 0')
        else:
           a = np.matrix('1. 0; 0 1.')

        self.A  = a
        self.dim = a.shape[1]

    def initialData(self, data):

        def u_a(x, t):   return       sin(2*pi*(x - t))
        def uT_a(x, t):  return -2*pi*cos(2*pi*(x - t))
        def uX_a(x, t):  return  2*pi*cos(2*pi*(x - t))

    def createData(self, grid_obj):

        def gE(x, t): return self.u_a (x, t)
        def gW(x, t): return self.u_a (x, t)
        def f(x, t):  return self.u_a(x, t)
