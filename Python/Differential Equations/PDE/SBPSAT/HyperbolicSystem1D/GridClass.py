class GridClass:
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
        print('plotting grid')

        plt.plot(np.zeros(self.Nx+1), self.grid)
        plt.show()
