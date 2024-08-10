class SchemeClass:
    def __init__(self, problem):
        self.A  = problem.A
        self.B  = problem.B

    def numericalBoundaryOperators(self, grid, problem):
        import numpy as np
        print('Creating numerical boundary operators')

        Ixy = np.eye((grid.Nx + 1) * (grid.Ny + 1))
        self.Ix = np.eye(grid.Nx + 1)
        self.Iy = np.eye(grid.Ny + 1)

        self.E0x = np.zeros((grid.Nx, grid.Nx))
        self.E0y = np.zeros((grid.Ny, grid.Ny))
        self.ENx = np.zeros((grid.Nx, grid.Nx))
        self.ENy = np.zeros((grid.Ny, grid.Ny))

        self.H_E_0_m_n = np.kron(np.kron(self.Ix, self.Iy),self.H_E_0_m_c)
        self.H_E_0_p_n = np.kron(np.kron(self.Ix, self.Iy),self.H_E_0_p_c)
        self.H_W_0_m_n = np.kron(np.kron(self.Ix, self.Iy),self.H_W_0_m_c)
        self.H_W_0_p_n = np.kron(np.kron(self.Ix, self.Iy),self.H_W_0_p_c)
        self.H_N_0_m_n = np.kron(np.kron(self.Ix, self.Iy),self.H_N_0_m_c)
        self.H_N_0_p_n = np.kron(np.kron(self.Ix, self.Iy),self.H_N_0_p_c)
        self.H_S_0_m_n = np.kron(np.kron(self.Ix, self.Iy),self.H_S_0_m_c)
        self.H_S_0_p_n = np.kron(np.kron(self.Ix, self.Iy),self.H_S_0_p_c)

        self.H_E_m_n = self.H_E_0_m_n
        self.H_W_p_n = self.H_W_0_p_n
        self.H_N_m_n = self.H_N_0_m_n
        self.H_S_p_n = self.H_S_0_p_n

        self.BC_E = self.H_E_m_n
        self.BC_W = self.H_W_p_n
        self.BC_N = self.H_N_m_n
        self.BC_S = self.H_S_p_n

        self.Lambda_E_m = -np.abs(np.linalg.solve(problem.A, np.eye(self.A.shape[1])))
        self.Lambda_E_m =  np.kron(np.kron(self.Ix, self.Iy), self.Lambda_E_m)
        self.Lambda_W_p =  np.linalg.solve(problem.A, np.eye(self.A.shape[1]))
        self.Lambda_W_p =  np.kron(np.kron(self.Ix, self.Iy), self.Lambda_W_p)

        self.Lambda_N_m = -np.abs(np.linalg.solve(problem.B, np.eye(self.B.shape[1])))
        self.Lambda_N_m =  np.kron(np.kron(self.Ix, self.Iy), self.Lambda_N_m)
        self.Lambda_S_p =  np.linalg.solve(problem.B, np.eye(self.B.shape[1]))
        self.Lambda_S_p =  np.kron(np.kron(self.Ix, self.Iy), self.Lambda_S_p)

    def ContinuousBoundaryOperators(self):
        import numpy as np
        print('Creating boundary operators')

        self.H_W_0_p_c = np.zeros((self.A.shape))
        self.H_W_0_m_c = np.zeros((self.A.shape))
        self.H_E_0_p_c = np.zeros((self.A.shape))
        self.H_E_0_m_c = np.zeros((self.A.shape))
        self.H_S_0_p_c = np.zeros((self.B.shape))
        self.H_S_0_m_c = np.zeros((self.B.shape))
        self.H_N_0_p_c = np.zeros((self.B.shape))
        self.H_N_0_m_c = np.zeros((self.B.shape))
        self.Lambda_W_p_c = np.zeros((self.A.shape))
        self.Lambda_E_m_c = np.zeros((self.A.shape))
        self.Lambda_S_p_c = np.zeros((self.B.shape))
        self.Lambda_N_m_c = np.zeros((self.B.shape))

        for k in range(0, self.A.shape[1]):
            if self.A[k,k] > 0:
                self.H_W_0_p_c[k, k]    = 1
                self.H_W_0_m_c[k, k]    = 0
                self.H_E_0_p_c[k, k]    = 1
                self.H_E_0_m_c[k, k]    = 0
                self.Lambda_W_p_c[k, k] = self.A[k, k]

            elif self.A[k,k] < 0:
                self.H_W_0_p_c[k, k]    = 0
                self.H_W_0_m_c[k, k]    = 1
                self.H_E_0_p_c[k, k]    = 0
                self.H_E_0_m_c[k, k]    = 1
                self.Lambda_E_m_c[k, k] = self.A[k, k]

            else:
                self.H_W_0_p_c[k, k] = 1
                self.H_W_0_m_c[k, k] = 1
                self.H_E_0_p_c[k, k] = 1
                self.H_E_0_m_c[k, k] = 1

        for k2 in range(0, self.B.shape[1]):
            if self.B[k2, k2] > 0:
                self.H_S_0_p_c[k2, k2]    = 1
                self.H_S_0_m_c[k2, k2]    = 0
                self.H_N_0_p_c[k2, k2]    = 1
                self.H_N_0_m_c[k2, k2]    = 0
                self.Lambda_S_p_c[k2, k2] = self.B[k2, k2]

            elif self.B[k2 ,k2] < 0:
                self.H_S_0_p_c[k2, k2]    = 0
                self.H_S_0_m_c[k2, k2]    = 1
                self.H_N_0_p_c[k2, k2]    = 0
                self.H_N_0_m_c[k2, k2]    = 1
                self.Lambda_N_m_c[k2, k2] = self.B[k2, k2]

            else:
                self.H_S_0_p_c[k2, k2] = 1
                self.H_S_0_m_c[k2, k2] = 1
                self.H_N_0_p_c[k2, k2] = 1
                self.H_N_0_m_c[k2, k2] = 1

    def createPenalties(self):
        import numpy as np
        print('Creating penalties')

        self.Sigma_E =  np.dot(self.H_E_m_n.T, self.Lambda_E_m)
        self.Sigma_W = -np.dot(self.H_W_p_n.T, self.Lambda_W_p)
        self.Sigma_N =  np.dot(self.H_N_m_n.T, self.Lambda_N_m)
        self.Sigma_S = -np.dot(self.H_S_p_n.T, self.Lambda_S_p)

    def CreateMatrices(self, grid, sbp, problem):
        print('Create matrices')
        import numpy as np

        self.Ix = np.eye(grid.Nx + 1)
        self.Iy = np.eye(grid.Ny + 1)

        E_0x         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_0y         = np.zeros((grid.Ny + 1, grid.Ny + 1))
        E_Nx         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_Ny         = np.zeros((grid.Ny + 1, grid.Ny + 1))
        E_0x[0, 0]   = 1
        E_0y[0, 0]   = 1
        E_Nx[-1, -1] = 1
        E_Ny[-1, -1] = 1

        self.Isys = np.eye(problem.dim)

        self.E_East  = np.kron(E_Nx, np.eye(grid.Ny+1))
        self.E_West  = np.kron(E_0x, np.eye(grid.Ny+1))
        self.E_North = np.kron(np.eye(grid.Nx+1), E_Ny)
        self.E_South = np.kron(np.eye(grid.Nx+1), E_0y)

        self.PinvENSigmaE = np.dot(np.dot(np.kron(np.kron(sbp.Pinv_x, self.Iy), self.Isys), self.Sigma_E), np.kron(self.E_East, self.Isys))
        self.PinvE0SigmaW = np.dot(np.dot(np.kron(np.kron(sbp.Pinv_x, self.Iy), self.Isys), self.Sigma_W), np.kron(self.E_West, self.Isys))
        self.PinvENSigmaN = np.dot(np.dot(np.kron(np.kron(self.Ix, sbp.Pinv_y), self.Isys), self.Sigma_N), np.kron(self.E_North, self.Isys))
        self.PinvE0SigmaS = np.dot(np.dot(np.kron(np.kron(self.Ix, sbp.Pinv_y), self.Isys), self.Sigma_S), np.kron(self.E_South, self.Isys))

        self.PinvE0SigmaWH0 = np.dot(self.PinvE0SigmaW, self.H_W_0_p_n)
        self.PinvENSigmaEH0 = np.dot(self.PinvENSigmaE, self.H_E_0_m_n)
        self.PinvE0SigmaSH0 = np.dot(self.PinvE0SigmaS, self.H_S_0_p_n)
        self.PinvENSigmaNH0 = np.dot(self.PinvENSigmaN, self.H_N_0_m_n)

        self.Amat = problem.A
        self.Bmat = problem.B

        Am = np.kron(np.kron(self.Ix, self.Iy), self.Amat)
        Bm = np.kron(np.kron(self.Ix, self.Iy), self.Bmat)

        Dx = np.kron(np.kron(sbp.D1x, self.Iy), self.Isys)
        Dy = np.kron(np.kron(self.Ix, sbp.D1y), self.Isys)

        self.A = -np.dot(Am, Dx) - np.dot(Bm, Dy) + np.dot(self.PinvE0SigmaW, self.BC_W) + np.dot(self.PinvENSigmaE, self.BC_E) + np.dot(self.PinvE0SigmaS, self.BC_S) + np.dot(self.PinvENSigmaN, self.BC_N)


        print(-np.dot(Am, Dx) - np.dot(Bm, Dy))
        print(np.dot(self.PinvE0SigmaW, self.BC_W))
        print(np.dot(self.PinvENSigmaE, self.BC_E))
        print(np.dot(self.PinvE0SigmaS, self.BC_S))
        print(np.dot(self.PinvENSigmaN, self.BC_N))
        print(self.A)

    def Scheme(self, grid, sbp, problem):
        print('Create scheme')
        import matplotlib.pyplot as plt
        import numpy as np

        X = grid.gridX
        Y = grid.gridY

        xlen     = grid.Nx + 1
        ylen     = grid.Ny + 1
        X = np.reshape(X, xlen*ylen)
        Y = np.reshape(Y, xlen*ylen)
        xzeroVec = grid.minx * np.ones((xlen, ylen))
        xzeroVec = np.reshape(xzeroVec, xlen*ylen)
        yzeroVec = grid.miny * np.ones((xlen, ylen))
        yzeroVec = np.reshape(yzeroVec, xlen*ylen)
        xoneVec  = grid.maxx * np.ones((xlen, ylen))
        xoneVec = np.reshape(xoneVec, xlen*ylen)
        yoneVec  = grid.maxy * np.ones((xlen, ylen))
        yoneVec = np.reshape(yoneVec, xlen*ylen)

        tone = np.ones((xlen, ylen))
        tone = np.reshape(tone, xlen*ylen)

        Amat  = problem.A
        Bmat  = problem.B

        def u_a(x, y, t):   return  [           np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t))          , np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t))]
        def uT_a(x, y, t):  return  [-2*np.pi * np.cos(2*np.pi*(x - t)) - 2*np.pi * np.cos(2*np.pi*(y - t)), -2*np.pi * np.cos(2*np.pi*(x - t)) - 2*np.pi * np.cos(2*np.pi*(y - t))]
        def uX_a(x, y, t):  return  [ 2*np.pi * np.cos(2*np.pi*(x - t)),  2*np.pi * np.cos(2*np.pi*(x - t))]
        def uY_a(x, y, t):  return  [ 2*np.pi * np.cos(2*np.pi*(y - t)),  2*np.pi * np.cos(2*np.pi*(y - t))]

        def gE(x, y, t): return u_a(x, y, t)
        def gW(x, y, t): return u_a(x, y, t)
        def gN(x, y, t): return u_a(x, y, t)
        def gS(x, y, t): return u_a(x, y, t)
        def f(x, y):     return u_a(x, y, 0*tone)

        def force(t): return np.reshape(np.transpose(uT_a(X, Y, t*tone)), xlen*ylen*problem.dim) + np.dot(np.kron(np.kron(self.Ix, self.Iy), Amat), np.reshape(np.transpose(uX_a(X, Y, t*tone)), xlen*ylen*problem.dim)) \
        + np.dot(np.kron(np.kron(self.Ix, self.Iy), Bmat), np.reshape(np.transpose(uY_a(X, Y, t*tone)), xlen*ylen*problem.dim))

        def cA(t, u): return np.squeeze(np.asarray(force(t*tone))) + np.squeeze(np.asarray(np.dot(self.A, u))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvE0SigmaWH0, np.reshape(np.transpose(gW(xzeroVec, Y, t*tone)), xlen*ylen*problem.dim)))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvENSigmaEH0, np.reshape(np.transpose(gE(xoneVec, Y,  t*tone)), xlen*ylen*problem.dim)))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvE0SigmaSH0, np.reshape(np.transpose(gS(X, yzeroVec, t*tone)), xlen*ylen*problem.dim)))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvENSigmaNH0, np.reshape(np.transpose(gN(X, yoneVec,  t*tone)), xlen*ylen*problem.dim))))

        u = np.zeros((xlen*ylen*problem.dim, grid.Nt + 1))

        F = np.squeeze(np.asarray(np.reshape(np.transpose(f(X, Y)), xlen*ylen*problem.dim)))
        u[:,0] = F
        up = u[:,0]
        dt = grid.dt

        for k in range(0, grid.Nt):

            tn = grid.t[k]
            k1 = cA(tn, up)
            k2 = cA(tn + dt/2, up + dt/2 * k1)
            k3 = np.reshape(cA(tn + dt/2, up + dt/2 * k2), xlen*ylen*problem.dim, 1)
            k4 = np.reshape(cA(tn + dt, up + dt * k3), xlen*ylen*problem.dim, 1)
            u[:, k+1] = np.reshape(up + dt/6*(k1 + 2*k2 + 2*k3 + k4), xlen*ylen*problem.dim, 1)
            up = u[:, k+1]

        self.u = u[:, -1]

        plt.plot(np.transpose(self.u))
        plt.show()
    def plotSolution(self, grid, problem):
        print('Plotting solution')
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import numpy as np

        X = grid.gridX
        Y = grid.gridY
        xlen = grid.Nx + 1
        ylen = grid.Ny + 1
        v = np.zeros((xlen, ylen, problem.dim))

        for k in range(0, problem.dim-1):
                print(len(self.u))
                print(k)
                print(problem.dim)
                print(len(self.u[k::problem.dim]))
                v[:,:,k] = np.reshape(self.u[k::problem.dim], [xlen, ylen])

        Z = v[:,:,0]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

        # Customize the z axis.
        ax.set_zlim(-1.01, 1.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
    def Error(self, sbp, grid, problem, t):
        print('Computing Error')
        import matplotlib.pyplot as plt
        import numpy as np

        X = grid.gridX
        Y = grid.gridY
        xlen = grid.Nx + 1
        ylen = grid.Ny + 1
        X = np.reshape(X, xlen*ylen)
        Y = np.reshape(Y, xlen*ylen)
        tone = 1*np.ones((xlen, ylen))
        tone = np.reshape(tone, xlen*ylen)

        def u_a(x, y, t):   return  [np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t)), np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t))]

        error = np.abs(np.reshape(np.transpose(u_a(X, Y, tone)), (xlen*ylen*problem.dim, 1)) - np.reshape(self.u, (xlen*ylen*problem.dim, 1)))

        plt.plot(np.reshape(np.transpose(u_a(X, Y, tone)), (xlen*ylen*problem.dim, 1)))
        plt.show()

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, np.kron(np.kron(sbp.Px, sbp.Py), self.Isys)), error))
        self.error_L2 = self.error_L2[0,0]

        print(self.error_L2)
