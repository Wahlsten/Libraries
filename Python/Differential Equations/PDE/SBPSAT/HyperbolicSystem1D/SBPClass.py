class SBPClass:
    def __init__(self, grid_obj):
        import numpy as np
        self.acc = 2

        N  = grid_obj.Nx
        dx = grid_obj.dx

        if self.acc == 2:
                d = 1/dx
                V = 1/dx*np.ones((1,N+1))
                V[0,0]  = V[0,0]*2
                V[0,-1] = V[0,-1]*2
                self.Pinv = np.diag(V[0])
                self.P    = np.diag(1./V[0])
                Z = np.zeros((1, N+1))
                O = np.ones((1, N))
                Q = 0.5*(np.diag(Z[0]) + np.diag(O[0], 1) - np.diag(O[0], -1))
                Q[0,0] = -0.5
                Q[-1,-1] = 0.5

                self.D1 = np.dot(self.Pinv, Q)
                self.D2 = np.dot(self.D1, self.D1)
                self.BS = self.D1
