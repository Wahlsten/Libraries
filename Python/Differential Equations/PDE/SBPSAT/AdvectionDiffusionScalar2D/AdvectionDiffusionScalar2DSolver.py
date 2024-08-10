class gridClass:
    'Common base class for all employees'
    def __init__(self, minv, maxv, N):
        import numpy as np

        self.minx = minv[0]
        self.miny = minv[1]
        self.maxx = maxv[0]
        self.maxy = maxv[1]
        self.Nx = N[0]
        self.Ny = N[1]
        self.Nt = N[2]
        self.dx = (maxv[0] - minv[0])/float(N[0])
        self.dy = (maxv[1] - minv[1])/float(N[1])
        self.dt = (maxv[2] - minv[2])/float(N[2])
        self.x = np.linspace(minv[0], maxv[0], self.Nx + 1, endpoint=True)
        self.y = np.linspace(minv[1], maxv[1], self.Ny + 1, endpoint=True)
        self.t = np.linspace(minv[2], maxv[2], self.Nt + 1, endpoint=True)
        self.gridy, self.gridx = np.meshgrid(self.y, self.x)

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print "plotting grid"

        '''plt.plot(np.zeros(self.Nx+1), self.grid)'''

class problemClass:
    def __init__(self, grid):
        import numpy as np

        def a_func (x, y):  return (10 + np.sin(x) * np.cos(y))
        def ax_func(x, y):  return (-np.cos(x) * np.cos(y))
        def b_func (x, y):  return (-2 + -np.cos(x) * np.sin(y))
        def by_func(x, y):  return (np.cos(x) * np.cos(y))
        def c_func (x, y):  return (0.01 + 0*x + 0*y)
        def d_func (x, y):  return (0.01 + 0*x + 0*y)
        def e_func (x, y):  return (0 + 0*x + 0*y)
        def f_func (x, y):  return (0 + 0*x + 0*y)

        A_m  = a_func(grid.gridx, grid.gridy)
        B_m  = b_func(grid.gridx, grid.gridy)
        C_m  = c_func(grid.gridx, grid.gridy)
        D_m  = d_func(grid.gridx, grid.gridy)
        E_m  = e_func(grid.gridx, grid.gridy)
        F_m  = f_func(grid.gridx, grid.gridy)
        Ax_m = ax_func(grid.gridx, grid.gridy)
        By_m = by_func(grid.gridx, grid.gridy)


        self.A  = np.diag(A_m.flatten(1))
        self.B  = np.diag(B_m.flatten(1))
        self.C  = np.diag(C_m.flatten(1))
        self.D  = np.diag(D_m.flatten(1))
        self.E  = np.diag(E_m.flatten(1))
        self.F  = np.diag(F_m.flatten(1))
        self.Ax = np.diag(Ax_m.flatten(1))
        self.By = np.diag(By_m.flatten(1))

    def initialData(self, data):

        def u_a(x, t):   return        sin(2*pi*(x - t))
        def uT_a(x, t):  return -2*pi* cos(2*pi*(x - t))
        def uX_a(x, t):  return  2*pi* cos(2*pi*(x - t))

    def createData(self, grid_obj):

        def gE(x, t): return self.u_a(x, t)
        def gW(x, t): return self.u_a(x, t)
        def f(x, t):  return self.u_a(x, t)

class SchemeClass:
    def __init__(self, problem):
        self.Am  = problem.A
        self.Bm  = problem.B
        self.Cm  = problem.C
        self.Dm  = problem.D
        self.Em  = problem.E
        self.Fm  = problem.F
        self.Axm  = problem.Ax
        self.Bym  = problem.By

    def numericalBoundaryOperators(self, grid, problem):
        import numpy as np
        print "Creating numerical boundary operators"

        Ix = np.eye((grid.Nx + 1))
        Iy = np.eye((grid.Ny + 1))
        Ixy = np.eye((grid.Nx + 1) * (grid.Ny + 1))

        Dx = np.kron(sbp.D1x, Iy)
        Dy = np.kron(Ix, sbp.D1y)

        self.E0x = np.zeros((grid.Nx, grid.Nx))
        self.E0y = np.zeros((grid.Ny, grid.Ny))
        self.ENx = np.zeros((grid.Nx, grid.Nx))
        self.ENy = np.zeros((grid.Ny, grid.Ny))

        self.H_E_0_m_n = self.H_E_0_m_c * Ixy
        self.H_E_0_p_n = self.H_E_0_p_c * Ixy
        self.H_W_0_m_n = self.H_W_0_m_c * Ixy
        self.H_W_0_p_n = self.H_W_0_p_c * Ixy
        self.H_N_0_m_n = self.H_N_0_m_c * Ixy
        self.H_N_0_p_n = self.H_N_0_p_c * Ixy
        self.H_S_0_m_n = self.H_S_0_m_c * Ixy
        self.H_S_0_p_n = self.H_S_0_p_c * Ixy

        self.H_E_Dx_p_n = np.dot(self.H_E_Dx_p_c, Dx)
        self.H_E_Dx_m_n = np.dot(self.H_E_Dx_m_c, Dx)
        self.H_W_Dx_p_n = np.dot(self.H_W_Dx_p_c, Dx)
        self.H_W_Dx_m_n = np.dot(self.H_W_Dx_m_c, Dx)
        self.H_N_Dx_p_n = np.dot(self.H_N_Dx_p_c, Dx)
        self.H_N_Dx_m_n = np.dot(self.H_N_Dx_m_c, Dx)
        self.H_S_Dx_p_n = np.dot(self.H_S_Dx_p_c, Dx)
        self.H_S_Dx_m_n = np.dot(self.H_S_Dx_m_c, Dx)

        self.H_E_Dy_p_n = np.dot(self.H_E_Dy_p_c, Dy)
        self.H_E_Dy_m_n = np.dot(self.H_E_Dy_m_c, Dy)
        self.H_W_Dy_p_n = np.dot(self.H_W_Dy_p_c, Dy)
        self.H_W_Dy_m_n = np.dot(self.H_W_Dy_m_c, Dy)
        self.H_N_Dy_p_n = np.dot(self.H_N_Dy_p_c, Dy)
        self.H_N_Dy_m_n = np.dot(self.H_N_Dy_m_c, Dy)
        self.H_S_Dy_p_n = np.dot(self.H_S_Dy_p_c, Dy)
        self.H_S_Dy_m_n = np.dot(self.H_S_Dy_m_c, Dy)

        self.H_W_p_n = (self.H_W_0_p_n + self.H_W_Dx_p_n + self.H_W_Dy_p_n)
        self.H_E_m_n = (self.H_E_0_m_n + self.H_E_Dx_m_n + self.H_E_Dy_m_n)
        self.H_S_p_n = (self.H_S_0_p_n + self.H_S_Dx_p_n + self.H_S_Dy_p_n)
        self.H_N_m_n = (self.H_N_0_m_n + self.H_N_Dx_m_n + self.H_N_Dy_m_n)

        self.H_E_m_n = self.H_E_0_m_n
        self.H_W_p_n = self.H_W_0_p_n
        self.H_N_m_n = self.H_N_0_m_n
        self.H_S_p_n = self.H_S_0_p_n

        self.BC_E = self.H_E_m_n
        self.BC_W = self.H_W_p_n
        self.BC_N = self.H_N_m_n
        self.BC_S = self.H_S_p_n

        self.Lambda_E_m = -np.abs(np.linalg.solve(problem.A, 1))
        self.Lambda_E_m =  np.diag(self.Lambda_E_m[:,0])
        self.Lambda_W_p =  np.linalg.solve(problem.A, 1)
        self.Lambda_W_p =  np.diag(self.Lambda_W_p[:,0])

        self.Lambda_N_m = -np.abs(np.linalg.solve(problem.B, 1))
        self.Lambda_N_m =  np.diag(self.Lambda_N_m[:,0])
        self.Lambda_S_p =  np.linalg.solve(problem.B, 1)
        self.Lambda_S_p =  np.diag(self.Lambda_S_p[:,0])

        print problem.A
    def ContinuousBoundaryOperators(self):
        import numpy as np
        print "Creating boundary operators"

        if self.Am[0, 0] > 0:

            self.H_W_0_p_c  = self.Am[0, 0]
            self.H_W_0_m_c  = 0

        elif self.Am[0, 0] < 0:

            self.H_W_0_p_c  = 0
            self.H_W_0_m_c  = self.Am[0, 0]

        else:

            self.H_W_0_p_c  = 1
            self.H_W_0_m_c  = 1

        if self.Am[-1, -1] > 0:

            self.H_E_0_p_c  = self.Am[-1, -1]
            self.H_E_0_m_c  = 0

        elif self.Am[-1, -1] < 0:

            self.H_E_0_p_c  = 0
            self.H_E_0_m_c  = self.Am[-1, -1]

        else:

            self.H_E_0_p_c  = 1
            self.H_E_0_m_c  = 1

        if self.Bm[0, 0] > 0:

            self.H_S_0_p_c  = self.Bm[0, 0]
            self.H_S_0_m_c  = 0

        elif self.Bm[0, 0] < 0:

            self.H_S_0_p_c  = 0
            self.H_S_0_m_c  = self.Bm[0, 0]

        else:

            self.H_S_0_p_c  = 1
            self.H_S_0_m_c  = 1

        if self.Bm[-1, -1] > 0:

            self.H_N_0_p_c  = self.Bm[-1, -1]
            self.H_N_0_m_c  = 0

        elif self.Bm[-1, -1] < 0:

            self.H_N_0_p_c  = 0
            self.H_N_0_m_c  = self.Bm[-1, -1]

        else:

            self.H_N_0_p_c  = 1
            self.H_N_0_m_c  = 1

        self.H_W_Dx_p_c = -self.Cm
        self.H_W_Dx_m_c = -self.Cm
        self.H_E_Dx_p_c = -self.Cm
        self.H_E_Dx_m_c = -self.Cm
        self.H_S_Dx_p_c = -self.Dm
        self.H_S_Dx_m_c = -self.Dm
        self.H_N_Dx_p_c = -self.Dm
        self.H_N_Dx_m_c = -self.Dm

        self.H_W_Dy_p_c = -self.Cm
        self.H_W_Dy_m_c = -self.Cm
        self.H_E_Dy_p_c = -self.Cm
        self.H_E_Dy_m_c = -self.Cm
        self.H_S_Dy_p_c = -self.Dm
        self.H_S_Dy_m_c = -self.Dm
        self.H_N_Dy_p_c = -self.Dm
        self.H_N_Dy_m_c = -self.Dm

    def createPenalties(self):
        import numpy as np
        print "Creating penalties"

        self.Sigma_E =  np.dot(self.H_E_m_n.T, self.Lambda_E_m)
        self.Sigma_W = -np.dot(self.H_W_p_n.T, self.Lambda_W_p)
        self.Sigma_N =  np.dot(self.H_N_m_n.T, self.Lambda_N_m)
        self.Sigma_S = -np.dot(self.H_S_p_n.T, self.Lambda_S_p)

        print self.Lambda_E_m
        print self.Lambda_W_p
        print self.Lambda_N_m
        print self.Lambda_S_p

    def CreateMatrices(self, grid, sbp, problem):
        print 'Create matrices'
        import numpy as np

        Ix = np.eye((grid.Nx + 1))
        Iy = np.eye((grid.Ny + 1))

        E_0x         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_Nx         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_0y         = np.zeros((grid.Ny + 1, grid.Ny + 1))
        E_Ny         = np.zeros((grid.Ny + 1, grid.Ny + 1))
        E_0y[ 0,  0] = 1
        E_Ny[-1, -1] = 1
        E_0x[ 0,  0] = 1
        E_Nx[-1, -1] = 1

        self.E_East  = np.kron(E_Nx, Iy)
        self.E_West  = np.kron(E_0x, Iy)
        self.E_North = np.kron(Ix, E_Ny)
        self.E_South = np.kron(Ix, E_0y)

        self.PinvENSigmaE = np.dot(np.dot(np.kron(sbp.Pinv_x, Iy), self.Sigma_E), self.E_East)
        self.PinvE0SigmaW = np.dot(np.dot(np.kron(sbp.Pinv_x, Iy), self.Sigma_W), self.E_West)
        self.PinvENSigmaN = np.dot(np.dot(np.kron(Ix, sbp.Pinv_y), self.Sigma_N), self.E_North)
        self.PinvE0SigmaS = np.dot(np.dot(np.kron(Ix, sbp.Pinv_y), self.Sigma_S), self.E_South)

        self.PinvE0SigmaWH0 = np.dot(self.PinvE0SigmaW, self.H_W_0_p_c)
        self.PinvENSigmaEH0 = np.dot(self.PinvENSigmaE, self.H_E_0_m_c)
        self.PinvE0SigmaSH0 = np.dot(self.PinvE0SigmaS, self.H_S_0_p_c)
        self.PinvENSigmaNH0 = np.dot(self.PinvENSigmaN, self.H_N_0_m_c)

        self.PinvE0SigmaWHDx = np.dot(self.PinvE0SigmaW, self.H_W_Dx_p_c)
        self.PinvENSigmaEHDx = np.dot(self.PinvENSigmaE, self.H_E_Dx_m_c)
        self.PinvE0SigmaSHDx = np.dot(self.PinvE0SigmaS, self.H_S_Dx_p_c)
        self.PinvENSigmaNHDx = np.dot(self.PinvENSigmaN, self.H_N_Dx_m_c)

        self.PinvE0SigmaWHDy = np.dot(self.PinvE0SigmaW, self.H_W_Dy_p_c)
        self.PinvENSigmaEHDy = np.dot(self.PinvENSigmaE, self.H_E_Dy_m_c)
        self.PinvE0SigmaSHDy = np.dot(self.PinvE0SigmaS, self.H_S_Dy_p_c)
        self.PinvENSigmaNHDy = np.dot(self.PinvENSigmaN, self.H_N_Dy_m_c)

        self.Amat  = problem.A
        self.Bmat  = problem.B
        self.Cmat  = problem.C
        self.Dmat  = problem.D
        self.Axmat = problem.Ax
        self.Bymat = problem.By

        Am   = self.Am
        Bm   = self.Bm
        Cm   = self.Cm
        Dm   = self.Dm
        Em   = self.Em
        Fm   = self.Fm
        Axm  = self.Axm
        Bym  = self.Bym

        Dx = np.kron(sbp.D1x, Iy)
        Dy = np.kron(Ix, sbp.D1y)

        self.A = - (np.dot(Am, Dx) + 0*np.dot(Dx, Am) - 0*Axm) \
                 - (np.dot(Bm, Dy) + 0*np.dot(Dy, Bm) - 0*Bym) \
                 + np.dot(np.dot(Cm, Dx), Dx) \
                 + np.dot(np.dot(Dm, Dy), Dy) \
                 + np.dot(np.dot(Em, Dx), Dy) \
                 + np.dot(np.dot(Fm, Dy), Dx) \
                 + np.dot(self.PinvE0SigmaW, self.BC_W) \
                 + np.dot(self.PinvENSigmaE, self.BC_E) \
                 + np.dot(self.PinvE0SigmaS, self.BC_S) \
                 + np.dot(self.PinvENSigmaN, self.BC_N)

    def Scheme(self, grid, sbp, problem):
        print 'Create scheme'
        import matplotlib.pyplot as plt
        import numpy as np

        X = grid.gridx
        Y = grid.gridy

        X = X.flatten()
        Y = Y.flatten()

        xlen     = grid.Nx + 1
        ylen     = grid.Ny + 1
        xzeroVec = grid.minx * np.ones((1, xlen))
        xzeroVec = xzeroVec[0, :]
        xoneVec  = grid.maxx * np.ones((1, xlen))
        xoneVec  = xoneVec[0, :]

        tone = np.ones((1, xlen*ylen))
        tone = tone[0, :]

        Amat  = problem.A
        Bmat  = problem.B
        Cmat  = problem.C
        Dmat  = problem.D
        Axmat = problem.Ax
        Bymat = problem.By

        def u_a(x, y, t):   return            np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t))
        def uT_a(x, y, t):  return -2*np.pi * np.cos(2*np.pi*(x - t)) - 2*np.pi * np.cos(2*np.pi*(y - t))
        def uX_a(x, y, t):  return  2*np.pi * np.cos(2*np.pi*(x - t)) + 0*y
        def uY_a(x, y, t):  return  2*np.pi * np.cos(2*np.pi*(y - t)) + 0*x
        def uXX_a(x, y, t): return -2*np.pi * 2*np.pi * np.sin(2*np.pi*(x - t)) + 0*y
        def uYY_a(x, y, t): return -2*np.pi * 2*np.pi * np.sin(2*np.pi*(y - t)) + 0*x

        def gE(x, y, t):   return u_a(x, y, t)
        def gW(x, y, t):   return u_a(x, y, t)
        def gN(x, y, t):   return u_a(x, y, t)
        def gS(x, y, t):   return u_a(x, y, t)

        def gEDx(x, y, t): return uX_a(x, y, t)
        def gWDx(x, y, t): return uX_a(x, y, t)
        def gNDx(x, y, t): return uX_a(x, y, t)
        def gSDx(x, y, t): return uX_a(x, y, t)

        def gEDy(x, y, t): return uY_a(x, y, t)
        def gWDy(x, y, t): return uY_a(x, y, t)
        def gNDy(x, y, t): return uY_a(x, y, t)
        def gSDy(x, y, t): return uY_a(x, y, t)

        def f(x, y):       return u_a(x, y, 0*tone)

        def force(t): return uT_a(X, Y, t*tone) \
        + np.dot(Amat, uX_a(X, Y,  t*tone)) \
        + np.dot(Bmat, uY_a(X, Y,  t*tone)) \
        - np.dot(Cmat, uXX_a(X, Y, t*tone)) \
        - np.dot(Dmat, uYY_a(X, Y, t*tone))

        def cA(t, u): return force(t*tone) + np.dot(self.A, u) \
        - np.dot(self.PinvE0SigmaWH0,  gW(  X, Y,  t*tone)) \
        - np.dot(self.PinvE0SigmaWHDx, gWDx(X, Y,  t*tone)) \
        - np.dot(self.PinvE0SigmaWHDy, gWDy(X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaEH0,  gE(  X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaEHDx, gEDx(X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaEHDy, gEDy(X, Y,  t*tone)) \
        - np.dot(self.PinvE0SigmaSH0,  gS(  X, Y,  t*tone)) \
        - np.dot(self.PinvE0SigmaSHDx, gSDx(X, Y,  t*tone)) \
        - np.dot(self.PinvE0SigmaSHDy, gSDy(X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaNH0,  gN(  X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaNHDx, gNDx(X, Y,  t*tone)) \
        - np.dot(self.PinvENSigmaNHDy, gNDy(X, Y,  t*tone))

        u = np.zeros((xlen * ylen, grid.Nt + 1))
        u[:, 0] = f(X, Y)
        up = u[:, 0]
        dt = grid.dt

        for k in range(0, grid.Nt):

            tn = grid.t[k]
            k1 = cA(tn, up)
            k2 = cA(tn + dt/2, up + dt/2 * k1)
            k3 = cA(tn + dt/2, up + dt/2 * k2)
            k4 = cA(tn + dt, up + dt * k3)
            u[:, k+1] = up + dt/6*(k1 + 2*k2 + 2*k3 + k4)
            up = u[:, k+1]

        self.u = u[:, -1]

    def Error(self, sbp, grid, problem):
        print 'Error'
        import numpy as np

        X = grid.gridx
        Y = grid.gridy
        X = X.flatten()
        Y = Y.flatten()

        def u_a(x, y, t): return np.sin(2*np.pi*(x - t)) + np.sin(2*np.pi*(y - t))

        error = np.abs(u_a(X, Y, grid.t[-1]) - self.u)

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, np.kron(sbp.Px, sbp.Py)), error))

        print self.error_L2

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

g = gridClass([0, 0, 0], [1, 1, 1], [20, 20, 4000])
g.plot_grid()

p = problemClass(g)
p.createData(g)

sbp     = SBPClass(g)
sbp.acc = 2

s = SchemeClass(p)
s.ContinuousBoundaryOperators()
s.numericalBoundaryOperators(g, p)
s.createPenalties()
s.CreateMatrices(g, sbp, p)
s.Scheme(g, sbp, p)
s.Error(sbp, g, p)
