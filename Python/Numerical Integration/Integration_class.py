class IntegrationClass:
    def __init__(self, minv, maxv, N, method, f, tol):
        import numpy as np
        self.minv = minv
        self.maxv = maxv
        self.N = N
        self.method = method
        self.func = f
        self.tol = tol

    def integrate(self):

        if self.method is 'MidpointRule':
            from MidpointRule import MidpointRule
            self.I = MidpointRule([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'TrapeziodalRule':
            from TrapeziodalRule import TrapeziodalRule
            self.I = TrapeziodalRule([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'SimpsonsRule':
            from SimpsonsRule import SimpsonsRule
            self.I = SimpsonsRule([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'BoolesRule':
            from BoolesRule import BoolesRule
            self.I = BoolesRule([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'GaussLegendre':
            from GaussLegendre import GaussLegendre
            self.I = GaussLegendre([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'GaussLobatto':
            from GaussLobatto import GaussLobatto
            self.I = GaussLobatto([self.minv, self.maxv], self.func, self.N)
        elif self.method is 'Romberg':
            from Romberg import Romberg
            self.I = Romberg([self.minv, self.maxv], self.func, self.N, self.tol)


#def f (x): return x*x*x*x

#Int = IntegrationClass(0, 1, 12, 'BoolesRule', f, 0.001)
#Int.integrate()
#print Int.I
