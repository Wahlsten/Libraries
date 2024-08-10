class gridClass:
    def __init__(self, minv, maxv, N, R):
        import numpy as np
        self.xmin = minv[0]
        self.xmax = maxv[0]
        self.Ns = N[0]
        self.Nt = N[1]
        self.dt = (maxv[1]-minv[1])/float(N[1])
        self.t = np.linspace(minv[1], maxv[1], N[1]+1, endpoint=True)
        self.s = np.sqrt(self.dt) * np.random.normal(0, 1, self.Ns * self.Nt)
        self.s = np.reshape(self.s, (self.Ns, self.Nt))
        self.R = R

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print("plotting grid")

        plt.plot(np.zeros(self.Nx+1), self.x)

class problemClass:
    def __init__(self, mu, sigma):
        import numpy as np

        self.mu = mu
        self.sigma = sigma

class SchemeClass:
    def __init__(self, problem, method):

        self.mu = problem.mu
        self.sigma = problem.sigma

        if method == 'EulerMaruyama':
            from EulerMaruyama import EulerMaruyama
            self.method = EulerMaruyama
            self.methodName = 'Euler-Maruyama'
        elif method == 'Milstein':
            from Milstein import Milstein
            self.method = Milstein
            self.methodName = 'Milstein'
        elif method == 'StochasticRungeKutta':
            from StochasticRungeKutta import StochasticRungeKutta
            self.method = StochasticRungeKutta
            self.methodName = 'StochasticRungeKutta'

    def Scheme(self, grid, problem):
        print('Create scheme')
        import matplotlib.pyplot as plt
        import numpy as np
        from math import pi
        from scipy import interpolate
        from matplotlib import animation

        solver = self.method

        ys = np.zeros((grid.Nt, 1))

        def afunc(x, t):   return 2
        def b_xfunc(x, t): return 1
        def bfunc(x, t):   return 2
        def f(t):          return 0

        Z       = grid.s
        a       = afunc
        b       = bfunc
        b_x     = b_xfunc
        R       = grid.R
        dt      = R * grid.dt * np.ones((1, grid.Ns))
        self.un = np.zeros((grid.Ns, int(grid.Nt/R)+1))
        self.un[:, 0] = f(0)

        for n in range(0,int(grid.Nt/R)):

            tp = grid.s[:, int(R*n)]
            Xp = self.un[:, n]
            Zp = np.cumsum(Z[:, int(R*n):int(R*n) + 1], axis = 1)

            self.un[:, n+1] = solver(self.un[:, n], a, b, b_x, dt, tp, Zp[:,-1])

        self.un_mean = np.mean(self.un, axis = 0)
        fig = plt.figure()
        ax = plt.axes(xlim=(grid.xmin,grid.xmax), ylim=(np.min(self.un), np.max(self.un)*1.1))
        plt.xlabel('t')
        plt.ylabel('u')

        for k in range(0, grid.Ns):

                plt.plot(grid.t[::int(R)], self.un[k, :])

        plt.show()

        fig2 = plt.figure()
        ax2 = plt.axes(xlim=(grid.xmin,grid.xmax), ylim=(np.min(self.un), np.max(self.un)*1.1))
        plt.xlabel('t')
        plt.ylabel('u')
        plt.plot(grid.t[::int(R)], self.un_mean)
        plt.show()

    def Error(self, grid, problem):
        print('Error')
        import numpy as np

        X = grid.x

        def u_a(x, t): return np.sin(2*np.pi*(x - t))

        error = np.abs(u_a(X, grid.t[-1]) - self.un)

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, sbp.P), error))

        print(self.error_L2)

g = gridClass([0, 0], [1, 1], [200, 200], 10.)

p = problemClass(0, 0)

s = SchemeClass(p, 'StochasticRungeKutta')
s.Scheme(g, p)
#s.Error(g, p)
