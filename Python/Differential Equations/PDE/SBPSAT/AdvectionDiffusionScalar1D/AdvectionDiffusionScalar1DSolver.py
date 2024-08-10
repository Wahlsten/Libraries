class gridClass:
    'Common base class for all employees'
    def __init__(self, minv, maxv, N):
        import numpy as np
        self.minx = minv[0]
        self.maxx = maxv[0]
        self.Nx = N[0]
        self.Nt = N[1]
        self.dx = (maxv[0]-minv[0])/float(N[0])
        self.dt = (maxv[1]-minv[1])/float(N[1])
        self.grid = np.linspace(minv[0], maxv[0], N[0] + 1, endpoint=True)
        self.t = np.linspace(minv[1], maxv[1], N[1]+1, endpoint=True)

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print "plotting grid"

        plt.plot(np.zeros(self.Nx+1), self.grid)

class problemClass:
    def __init__(self, grid):
        import numpy as np

        def a_func (x):  return (2 + 1*x)
        def ax_func(x):  return (1 + 0*x)

        self.A  = np.diag(a_func(grid.grid))
        self.Ax = np.diag(ax_func(grid.grid))

    def initialData(self, data):

        def u_a(x, t):   return           sin(2*pi*(x - t))
        def uT_a(x, t):  return -2*pi*    cos(2*pi*(x - t))
        def uX_a(x, t):  return  2*pi*    cos(2*pi*(x - t))

    def createData(self, grid_obj):

        def gE(x, t): return self.u_a (x, t)
        def gW(x, t): return self.u_a (x, t)
        def f(x, t): return self.u_a(x, t)

class SchemeClass:
    def __init__(self, problem):
        self.A  = problem.A

    def numericalBoundaryOperators(self, grid, problem):
        import numpy as np
        print "Creating numerical boundary operators"

        self.E0 = np.zeros((grid.Nx, grid.Nx))
        self.EN = np.zeros((grid.Nx, grid.Nx))

        self.H_E_0_m_n  = self.H_E_0_m_c * np.eye(grid.Nx + 1)
        self.H_E_0_p_n  = self.H_E_0_p_c * np.eye(grid.Nx + 1)
        self.H_W_0_m_n  = self.H_W_0_m_c * np.eye(grid.Nx + 1)
        self.H_W_0_p_n  = self.H_W_0_p_c * np.eye(grid.Nx + 1)

        self.H_W_p_n = self.H_W_0_p_n
        self.H_E_m_n = self.H_E_0_m_n

        self.BC_W = self.H_W_p_n
        self.BC_E = self.H_E_m_n

        self.Lambda_E_m = -np.linalg.solve(problem.A, 1)
        self.Lambda_E_m =  np.diag(self.Lambda_E_m[:,0])

        self.Lambda_W_p = np.linalg.solve(problem.A, 1)
        self.Lambda_W_p = np.diag(self.Lambda_W_p[:,0])

    def ContinuousBoundaryOperators(self):
        import numpy as np
        print "Creating boundary operators"

        if self.A[0, 0] > 0:

            self.H_W_0_p_c  = self.A[0, 0]
            self.H_W_0_m_c  = 0

        elif self.A[0, 0] < 0:

            self.H_W_0_p_c  = 0
            self.H_W_0_m_c  = self.A[0, 0]

        else:

            self.H_W_0_p_c  = 1
            self.H_W_0_m_c  = 1

        if self.A[-1, -1] > 0:

            self.H_E_0_p_c  = self.A[-1, -1]
            self.H_E_0_m_c  = 0

        elif self.A[-1, -1] < 0:

            self.H_E_0_p_c  = 0
            self.H_E_0_m_c  = self.A[-1, -1]

        else:

            self.H_E_0_p_c  = 1
            self.H_E_0_m_c  = 1

    def createPenalties(self):
        import numpy as np
        print "Creating penalties"

        self.Sigma_E =  np.dot(self.H_E_m_n.T, self.Lambda_E_m)
        self.Sigma_W = -np.dot(self.H_W_p_n.T, self.Lambda_W_p)

    def CreateMatrices(self, grid, sbp, problem):
        print 'Create matrices'
        import numpy as np

        E_0         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_N         = np.zeros((grid.Nx + 1, grid.Nx + 1))
        E_0[0, 0]   = 1
        E_N[-1, -1] = 1

        self.E_East = E_N
        self.E_West = E_0

        self.PinvENSigmaE = np.dot(np.dot(sbp.Pinv, self.Sigma_E), self.E_East)
        self.PinvE0SigmaW = np.dot(np.dot(sbp.Pinv, self.Sigma_W), self.E_West)

        self.PinvE0SigmaWH0  = np.dot(self.PinvE0SigmaW, self.H_W_0_p_c)
        self.PinvENSigmaEH0  = np.dot(self.PinvENSigmaE, self.H_E_0_m_c)

        self.Amat  = problem.A
        self.Axmat = problem.Ax

        Am   = self.Amat
        Axm  = self.Axmat

        Dx = sbp.D1

        self.A = -(1/2.)*(np.dot(Am, Dx) + np.dot(Dx, Am) - Axm) + np.dot(self.PinvE0SigmaW, self.BC_W) + np.dot(self.PinvENSigmaE, self.BC_E)

    def Scheme(self, grid, sbp, problem):
        print 'Create scheme'
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
        Axmat = problem.Ax

        '''def u_a(x, t):   return  1 + 2*x + 0*t**2
        def uT_a(x, t):  return  0 + 0*x + 0*t
        def uX_a(x, t):  return  2 + 0*x + 0*t
        def uXX_a(x, t): return  0 + 0*x + 0*t'''

        def u_a(x, t):   return             np.sin(2*np.pi*(x - t))
        def uT_a(x, t):  return  -2*np.pi * np.cos(2*np.pi*(x - t))
        def uX_a(x, t):  return   2*np.pi * np.cos(2*np.pi*(x - t))

        def gE(x, t):  return u_a(x, t)
        def gW(x, t):  return u_a(x, t)
        def f(x):      return u_a(x, 0*tone)

        def force(t): return uT_a(X, t*tone) + np.dot(Amat, uX_a(X, t*tone))

        def cA(t, u): return force(t*tone) + np.dot(self.A, u) - np.dot(self.PinvE0SigmaWH0, gW(xzeroVec, t*tone)) - np.dot(self.PinvENSigmaEH0, gE(xoneVec,  t*tone))

        u = np.zeros((xlen, grid.Nt+1))
        u[:,0] = f(X)
        up = u[:,0]
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

        plt.plot(grid.grid, self.u)
        plt.show()

    def Error(self, sbp, grid, problem):
        print 'Error'
        import numpy as np

        X = grid.grid

        def u_a(x, t): return np.sin(2*np.pi*(x - t))
        '''def u_a(x, t):   return  1 + 2*x + 0*t**2'''

        error = np.abs(u_a(X, grid.t[-1]) - self.u)

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, sbp.P), error))

        print self.error_L2

class SBPClass:
    def __init__(self, grid_obj):
        import numpy as np
        self.acc = 2

        N  = grid_obj.Nx
        dx = grid_obj.dx

        if self.acc == 2:
                d = 1/dx
                V = 1/dx*np.ones((1,N+1))
                V[0,0] = V[0,0]*2
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

g = gridClass([0, 0], [1, 1], [160, 4000])
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
