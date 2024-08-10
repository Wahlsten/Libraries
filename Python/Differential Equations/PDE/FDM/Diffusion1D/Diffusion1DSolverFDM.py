class gridClass:
    'Common base class for all employees'
    def __init__(self, minv, maxv, N):
        import numpy as np
        self.xmin = minv[0]
        self.xmax = maxv[0]
        self.Nx = N[0]
        self.Nt = N[1]
        self.dx = (maxv[0]-minv[0])/float(N[0])
        self.dt = (maxv[1]-minv[1])/float(N[1])
        self.x = np.linspace(minv[0], maxv[0], N[0] + 1, endpoint=True)
        self.t = np.linspace(minv[1], maxv[1], N[1]+1, endpoint=True)

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print "plotting grid"

        plt.plot(np.zeros(self.Nx+1), self.x)

class problemClass:
    def __init__(self, b):
        import numpy as np

        self.B = b

class SchemeClass:
    def __init__(self, problem, method):

        self.B  = problem.B

        if method is 'CrankNicolson':
            from CrankNicolson import CrankNicolson
            self.method = CrankNicolson
            self.methodname = 'Crank Nicolson'
            self.methodSolver = 'Implicit'
        elif method is 'ForwardTimeCentralSpace':
            from ForwardTimeCentralSpace import ForwardTimeCentralSpace
            self.method = ForwardTimeCentralSpace
            self.methodname = 'Forward time central space'
            self.methodSolver = 'Explicit'
        elif method is 'BackwardTimeCentralSpace':
            from BackwardTimeCentralSpace import BackwardTimeCentralSpace
            self.method = BackwardTimeCentralSpace
            self.methodname = 'Backward time central space'
            self.methodSolver = 'Implicit'
        elif method is 'DuFortFrankel':
            from DuFortFrankel import DuFortFrankel
            self.method = DuFortFrankel
            self.methodname = 'Du Fort Frankel'
            self.methodSolver = 'Multistep'

    def Scheme(self, grid, problem):
        print 'Create scheme'
        import matplotlib.pyplot as plt
        import numpy as np
        from math import pi
        from scipy import interpolate
        from matplotlib import animation

        solver = self.method
        B = problem.B
        def u_a(x, t): return np.sin(2*pi*x)*np.exp(-B*2*pi*2*pi*t)
        def g(x,t):      return u_a(x, t)*0
        def f(x):      return u_a(x, 0)


        self.un = np.zeros((grid.Nx + 1, grid.Nt + 1))
        self.un[:, 0] = u_a(grid.x, 0)

        C = grid.dt * B / (grid.dx * grid.dx)

        if self.methodSolver is 'Explicit':

            for n in range(0, grid.Nt):

                self.un[1:-1, n+1] = solver(self.un[:, n], C)

        elif self.methodSolver is 'Implicit':

                [D, H] = solver(grid.Nx, grid.Nt, grid.x, grid.t, B, grid.dt, grid.dx, f, g)

                self.un = np.linalg.solve(D,H)
                self.un = self.un.reshape(grid.Nt+1, grid.Nx+1)
                self.un = np.transpose(self.un)

        elif self.methodSolver is 'Multistep':

                self.un[1:-1, 1] = solver(self.un[:, 0], self.un[:, 0], C, 1)

                for n in range(1, grid.Nt):

                    self.un[1:-1, n+1] = solver(self.un[:, n], self.un[:, n-1], C, 0)
        plt.show()
        fig = plt.figure()
        ax = plt.axes(xlim=(grid.xmin,grid.xmax), ylim=(np.min(self.un), np.max(self.un)*1.1))

        lines=[]     # list for plot lines for solvers and analytical solutions
        legends=[]   # list for legends for solvers and analytical solutions


        line, = ax.plot([], [])
        lines.append(line)
        legends.append(self.methodname)

        line, = ax.plot([], []) #add extra plot line for analytical solution
        lines.append(line)
        legends.append('Analytical')

        plt.xlabel('x')
        plt.ylabel('u')
        plt.legend(legends, loc=1, frameon=False)

        # initialization function: plot the background of each frame
        def init():
            for line in lines:
                line.set_data([], [])
            return lines,

        # animation function.  This is called sequentially
        def animate(i):
            for k, line in enumerate(lines):
                if (k==0):
                    line.set_data(grid.x, self.un[:, i])
                else:
                    line.set_data(grid.x, u_a(grid.x, i*grid.dt))
            return lines,

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=grid.Nt, interval=100, blit=False)

        plt.show()
    def Error(self, grid, problem):
        print 'Error'
        import numpy as np

        X = grid.x

        def u_a(x, t): return np.sin(2*np.pi*(x - t))

        error = np.abs(u_a(X, grid.t[-1]) - self.un)

        self.error_L2 = np.sqrt(np.dot(np.dot(error.T, sbp.P), error))

        print self.error_L2

g = gridClass([0, 0], [1, 1], [50, 60])

p = problemClass(0.1)

s = SchemeClass(p, 'BackwardTimeCentralSpace')
s.Scheme(g, p)
#s.Error(g, p)
