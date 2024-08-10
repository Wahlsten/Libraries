class SchemeClass:
    def __init__(self, problem):
        self.A  = problem.A

    def numericalBoundaryOperators(self, grid, problem):
        import numpy as np
        print('Creating numerical boundary operators')

        self.Ix = np.eye(grid.Nx + 1)

        self.E0 = np.zeros((grid.Nx+1, grid.Nx+1))
        self.EN = np.zeros((grid.Nx+1, grid.Nx+1))

        self.H_E_0_m_n = np.kron(np.eye(grid.Nx + 1), self.H_E_0_m_c)
        self.H_E_0_p_n = np.kron(np.eye(grid.Nx + 1), self.H_E_0_p_c)
        self.H_W_0_m_n = np.kron(np.eye(grid.Nx + 1), self.H_W_0_m_c)
        self.H_W_0_p_n = np.kron(np.eye(grid.Nx + 1), self.H_W_0_p_c)

        self.H_W_p_n = self.H_W_0_p_n
        self.H_E_m_n = self.H_E_0_m_n

        self.BC_W = self.H_W_p_n
        self.BC_E = self.H_E_m_n

        '''
        self.Lambda_E_m = -np.linalg.solve(problem.A, 1)
        self.Lambda_E_m =  np.diag(self.Lambda_E_m[:,0])
        self.Lambda_E_m =  np.kron(self.Ix, self.Lambda_E_m)

        self.Lambda_W_p = np.linalg.solve(problem.A, 1)
        self.Lambda_W_p = np.diag(self.Lambda_W_p[:,0])
        self.Lambda_W_p = np.kron(self.Ix, self.Lambda_W_p)
        '''

        self.Lambda_E_m = np.kron(self.Ix, self.Lambda_E_m_c)
        self.Lambda_W_p = np.kron(self.Ix, self.Lambda_W_p_c)

        self.H_W_0_p_c = np.kron(self.Ix, self.H_W_0_p_c)
        self.H_E_0_m_c = np.kron(self.Ix, self.H_E_0_m_c)

    def ContinuousBoundaryOperators(self):
        import numpy as np
        print('Creating boundary operators')

        self.H_W_0_p_c = np.zeros((self.A.shape))
        self.H_W_0_m_c = np.zeros((self.A.shape))
        self.H_E_0_p_c = np.zeros((self.A.shape))
        self.H_E_0_m_c = np.zeros((self.A.shape))
        self.Lambda_W_p_c = np.zeros((self.A.shape))
        self.Lambda_E_m_c = np.zeros((self.A.shape))

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


    def createPenalties(self):
        import numpy as np
        print('Creating penalties')

        self.Sigma_E =  np.dot(self.H_E_m_n.T, self.Lambda_E_m)
        self.Sigma_W = -np.dot(self.H_W_p_n.T, self.Lambda_W_p)

    def CreateMatrices(self, grid, sbp, problem):
        print('Create matrices')
        import numpy as np

        E_0         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_N         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_0[0, 0]   = 1
        E_N[-1, -1] = 1

        self.Isys = np.eye(problem.dim)

        self.E_East = E_N
        self.E_West = E_0

        self.PinvENSigmaE = np.dot(np.dot(np.kron(sbp.Pinv, self.Isys), self.Sigma_E), np.kron(self.E_East, self.Isys))
        self.PinvE0SigmaW = np.dot(np.dot(np.kron(sbp.Pinv, self.Isys), self.Sigma_W), np.kron(self.E_West, self.Isys))

        self.PinvE0SigmaWH0 = np.dot(self.PinvE0SigmaW, self.H_W_0_p_c)
        self.PinvENSigmaEH0 = np.dot(self.PinvENSigmaE, self.H_E_0_m_c)

        self.Amat = problem.A

        Am = np.kron(self.Ix, self.Amat)

        Dx = np.kron(sbp.D1, self.Isys)
        self.A = -np.dot(Am, Dx) + np.dot(self.PinvE0SigmaW, self.BC_W) + np.dot(self.PinvENSigmaE, self.BC_E)

    def Scheme(self, grid, sbp, problem):
        print('Create scheme')
        import matplotlib.pyplot as plt
        import numpy as np

        X = grid.grid

        xlen     = grid.Nx + 1
        xzeroVec = grid.minx * np.ones((1, xlen))
        xzeroVec = xzeroVec[0,:]
        xoneVec  = grid.maxx * np.ones((1, xlen))
        xoneVec  = xoneVec[0,:]

        tone = np.ones((1, xlen))
        tone = tone[0,:]

        Amat  = problem.A

        def u_a(x, t):   return  [           np.sin(2*np.pi*(x - t)),            np.sin(2*np.pi*(x - t))]
        def uT_a(x, t):  return  [-2*np.pi * np.cos(2*np.pi*(x - t)), -2*np.pi * np.cos(2*np.pi*(x - t))]
        def uX_a(x, t):  return  [ 2*np.pi * np.cos(2*np.pi*(x - t)),  2*np.pi * np.cos(2*np.pi*(x - t))]

        def gE(x, t):  return u_a(x, t)
        def gW(x, t):  return u_a(x, t)
        def f(x):      return u_a(x, 0*tone)

        def force(t): return np.reshape(np.transpose(uT_a(X, t*tone)), xlen*problem.dim) + np.dot(np.kron(self.Ix, Amat), np.reshape(np.transpose(uX_a(X, t*tone)), xlen*problem.dim))

        def cA(t, u): return np.squeeze(np.asarray(force(t*tone))) + np.squeeze(np.asarray(np.dot(self.A, u))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvE0SigmaWH0, np.reshape(np.transpose(gW(xzeroVec, t*tone)), xlen*problem.dim)))) \
        -  np.squeeze(np.asarray(np.dot(self.PinvENSigmaEH0, np.reshape(np.transpose(gE(xoneVec,  t*tone)), xlen*problem.dim))))

        u = np.zeros((xlen*problem.dim, grid.Nt + 1))
        F = np.squeeze(np.asarray(np.reshape(np.transpose(f(X)), xlen*problem.dim)))
        u[:,0] = F
        up = u[:,0]
        dt = grid.dt

        for k in range(0, grid.Nt):

            tn = grid.t[k]
            k1 = cA(tn, up)
            k2 = cA(tn + dt/2, up + dt/2 * k1)
            k3 = np.reshape(cA(tn + dt/2, up + dt/2 * k2), xlen*problem.dim, 1)
            k4 = np.reshape(cA(tn + dt, up + dt * k3), xlen*problem.dim, 1)
            u[:, k+1] = np.reshape(up + dt/6*(k1 + 2*k2 + 2*k3 + k4), xlen*problem.dim, 1)
            up = u[:, k+1]

        self.u = u[:, -1]
        self.ut = u

    def plotSolution(self, grid, problem):
        print('Plotting solution')
        import numpy as np
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as p3
        import matplotlib.animation as animation

        X = grid.grid
        T = grid.t
        xlen = grid.Nx + 1
        tlen = grid.Nt + 1
        v = np.zeros((xlen, tlen, problem.dim))

        for k in range(0, problem.dim):

                v[:, :, k] = np.reshape(self.ut[k::problem.dim, :], [xlen, tlen])

        Y = v[:,:,0]
        Z = v[:,:,1]

        def update_lines(num, dataLines, lines):
            for line, data in zip(lines, dataLines):
                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, :num])
                line.set_3d_properties(data[2, :num])
            return lines

        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        # Fifty lines of random 3-D lines
        # data = [np.array([X, Y[:, 0], 0*Z[:,0]]), np.array([X, 0*Y[:, 0], Z[:,0]])]
        jump = 2
        data = [np.array([T[::jump], Y[0, ::jump], 0*Z[0,::jump]]), np.array([T[::jump], 0*Y[0, ::jump], Z[0, ::jump]])]

        lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

        # Setting the axes properties
        ax.set_xlim3d([0.0, T[-1]])
        ax.set_xlabel('T')

        ax.set_ylim3d([-1.0, 1.0])
        ax.set_ylabel('Y')

        ax.set_zlim3d([-1.0, 1.0])
        ax.set_zlabel('Z')

        ax.set_title('3D Test')

        # Creating the Animation object
        line_ani = animation.FuncAnimation(fig, update_lines, grid.Nt+1, fargs=(data, lines),
                                           interval=1, blit=False)

        plt.show()

    def Error(self, sbp, grid, problem, t):
        print('Computing Error')
        import numpy as np

        X = grid.grid
        xlen = grid.Nx + 1
        tone = np.ones((1, xlen))
        tone = tone[0,:]

        def u_a(x, t): return [np.sin(2*np.pi*(x - t)), np.sin(2*np.pi*(x - t))]

        error = np.abs(np.reshape(np.transpose(u_a(X, t*tone)), (xlen*problem.dim, 1)) - np.reshape(self.u, (xlen*problem.dim, 1)))

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, np.kron(sbp.P, self.Isys)), error))
        self.error_L2 = self.error_L2[0,0]

        print(self.error_L2)
