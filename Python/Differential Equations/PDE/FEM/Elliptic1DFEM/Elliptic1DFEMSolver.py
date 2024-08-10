class gridClass:
    'Common base class for all employees'
    def __init__(self, minv, maxv, N):
        import numpy as np
        self.minx = minv[0]
        self.maxx = maxv[0]
        self.Nx = N[0]
        self.Nt = N[1]
        '''self.grid = np.array([0., 0.1, 0.3, 0.33, 0.5, 0.75, 1.0])'''
        self.grid = np.linspace(minv[0], maxv[0], N[0] + 1, endpoint=True)
        self.dx = self.grid[1:] - self.grid[0:-1]
        self.dt = (maxv[1]-minv[1])/float(N[1])

        self.t = np.linspace(minv[1], maxv[1], N[1]+1, endpoint=True)

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print "plotting grid"

        plt.plot(np.zeros(self.Nx+1), self.grid)

class problemClass:
    def __init__(self, grid):
        import numpy as np

        def a_func(x): return (3 + 0*x)
        def b_func(x): return (1 + 0*x)

        self.a = a_func
        self.b = b_func
        self.A = np.diag(a_func(grid.grid))
        self.B = np.diag(b_func(grid.grid))

    def initialData(self, data):

        def u_a(x, t):   return        sin(2*pi*(x - t))
        def uT_a(x, t):  return -2*pi* cos(2*pi*(x - t))
        def uX_a(x, t):  return  2*pi* cos(2*pi*(x - t))

    def createData(self, grid_obj):

        def gE(x, t): return self.u_a (x, t)
        def gW(x, t): return self.u_a (x, t)
        def f(x, t):  return self.u_a(x, t)

class SchemeClass:
    def __init__(self, problem):
        self.A  = problem.A

    def Scheme(self, grid, problem):
        print 'Create scheme'
        import matplotlib.pyplot as plt
        import numpy as np

        def hat1(x,x1,x2): return (x - x1) / (x2 - x1)
        def hat2(x,x1,x2): return (x2 - x) / (x2 - x1)

        def u_a(x):    return np.sin(np.pi*x)
        def u_a_x(x):  return np.pi * np.cos(np.pi*x)
        def u_a_xx(x): return -np.pi*np.pi * np.sin(np.pi*x)

        def f(x):      return -problem.a(x)*u_a_xx(x) + problem.b(x)*u_a(x)

        x = grid.grid
        dx = grid.dx
        Nx = grid.Nx
        A = np.zeros((Nx+1, Nx+1))
        M = np.zeros((Nx+1, Nx+1))
        F = np.zeros((Nx+1, 1))
        A[0, 0]   = 1
        F[0]      = u_a(x[0])
        A[-1, -1] = 1
        F[-1]     = u_a(x[-1])
        A[1, 1]   = 1/dx[0]
        M[1, 1]   = dx[0]/3
        F[1]      = self.int_hat1_f(x[0],x[1], f)

        for i in range(1, Nx-1):
            A[i, i]     = A[i, i] + 1/dx[i]
            A[i, i+1]   = A[i, i+1] - 1/dx[i]
            A[i+1, i]   = A[i+1, i] - 1/dx[i]
            A[i+1, i+1] = A[i+1, i+1] + 1/dx[i]
            M[i, i]     = M[i, i] + dx[i]/3
            M[i, i+1]   = M[i, i+1] + dx[i]/6
            M[i+1, i]   = M[i+1, i] + dx[i]/6
            M[i+1, i+1] = M[i+1, i+1] + dx[i]/3
            F[i]        = F[i]   + self.int_hat2_f(x[i],x[i+1], f)
            F[i+1]      = F[i+1] + self.int_hat1_f(x[i],x[i+1], f)

        A[Nx-1, Nx-1] = A[Nx-1, Nx-1] + 1/dx[-1]
        M[Nx-1, Nx-1] = M[Nx-1, Nx-1] + dx[-1]/3

        F[Nx-1] = F[Nx-1] + self.int_hat2_f(x[Nx-1],x[Nx], f)
        Amat = np.dot(problem.A, A) + np.dot(problem.B, M)
        self.u = np.squeeze(np.asarray(np.linalg.solve(Amat, F)))

        plt.plot(grid.grid, self.u)
        plt.axis([0,1,0, np.max(self.u)])
        plt.show()

    def int_hat1_f(self, x1, x2, f):
        import numpy as np
        def hat1(x,x1,x2): return (x - x1) / (x2 - x1)
        def func(x,x1,x2): return hat1(x, x1, x2) * f(x)

        Int = self.SimpsonsRule(x1, x2, 10, func)

        return Int

    def int_hat2_f(self, x1, x2, f):
        import numpy as np
        def hat2(x,x1,x2): return (x2 - x) / (x2 - x1)
        def func(x,x1,x2): return hat2(x, x1, x2) * f(x)

        Int = self.SimpsonsRule(x1, x2, 10, func)

        return Int

    def fem_sol(self, x, x2):
        import numpy as np
        def hat1(x,x1,x2): return (x - x1) / (x2 - x1)
        def hat2(x,x1,x2): return (x2 - x) / (x2 - x1)

        for k in range(0,len(x)):
            if x[k] > x2:
                return hat2(x2,x[k-1],x[k])*self.u[k-1] + hat1(x2,x[k-1],x[k])*self.u[k]

        return hat2(x2,x[-2],x[-1])*self.u[-2] + hat1(x2,x[-2],x[-1])*self.u[-1]

    def SimpsonsRule(self, minv, maxv, N, func):
        import numpy as np

        grid = np.linspace(minv, maxv, N + 1, endpoint=True)
        dx = grid[1] - grid[0]

        I = 2*dx/6. * (func(grid[0], minv, maxv) + 4*np.sum(func(grid[1:-1:2], minv, maxv)) + 2*np.sum(func(grid[2:-1:2], minv, maxv)) + func(grid[-1], minv, maxv))

        return I

    def Error(self, grid, problem):
        print 'Error'
        import matplotlib.pyplot as plt
        import numpy as np

        N = 20
        xv = grid.grid
        xn = np.linspace(grid.minx, grid.maxx, N + 1, endpoint=True)
        self.un = np.squeeze(np.asarray(np.zeros((1, N+1))))

        for k in range(0, N):
            self.un[k] = self.fem_sol(xv, xn[k])


        def u_a(x):    return np.sin(np.pi*x)

        error = np.abs(u_a(xn) - self.un)
        self.error_L2 = np.sqrt(np.dot(error.T, error))

        print self.error_L2

        plt.plot(xn, u_a(xn), xn, self.un)
        plt.axis([0,1,0,np.max([np.max(u_a(xn)), np.max(self.un)])])
        plt.show()

        plt.plot(xn, error)
        plt.axis([grid.minx, grid.maxx,0,np.max(error)])
        plt.show()

g = gridClass([0, 0], [1, 1], [20, 4000])
g.plot_grid()

p = problemClass(g)
p.createData(g)

s = SchemeClass(p)
s.Scheme(g, p)
s.Error(g, p)
