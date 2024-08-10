class SBPClass:
    def __init__(self, grid_obj):
        import numpy as np
        self.acc = 2

        Nx = grid_obj.Nx
        Ny = grid_obj.Ny
        dx = grid_obj.dx
        dy = grid_obj.dy

        if self.acc == 2:
                d = 1/dx
                V = 1/dx*np.ones((1,Nx + 1))
                V[0, 0] = V[0,0]*2
                V[0,-1] = V[0,-1]*2
                self.Pinv_x = np.diag(V[0])
                self.Px     = np.diag(1./V[0])
                Z = np.zeros((1, Nx + 1))
                O = np.ones((1, Nx))
                Q = 0.5*(np.diag(Z[0]) + np.diag(O[0], 1) - np.diag(O[0], -1))
                Q[0,  0] = -0.5
                Q[-1,-1] =  0.5

                self.D1x = np.dot(self.Pinv_x, Q)
                self.D2x = np.dot(self.D1x, self.D1x)
                self.BSx = self.D1x

                d = 1/dy
                V = 1/dy * np.ones((1,Ny + 1))
                V[0, 0] = V[0,0]*2
                V[0,-1] = V[0,-1]*2
                self.Pinv_y = np.diag(V[0])
                self.Py     = np.diag(1./V[0])
                Z = np.zeros((1, Ny + 1))
                O = np.ones((1, Ny))
                Q = 0.5*(np.diag(Z[0]) + np.diag(O[0], 1) - np.diag(O[0], -1))
                Q[0,  0] = -0.5
                Q[-1,-1] =  0.5

                self.D1y = np.dot(self.Pinv_y, Q)
                self.D2y = np.dot(self.D1y, self.D1y)
                self.BSy = self.D1y
