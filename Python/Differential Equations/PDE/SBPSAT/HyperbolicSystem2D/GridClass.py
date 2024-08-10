class GridClass:
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
        self.dx = (maxv[0]-minv[0])/float(N[0])
        self.dy = (maxv[1]-minv[1])/float(N[1])
        self.dt = (maxv[2]-minv[2])/float(N[2])
        self.gridx = np.linspace(minv[0], maxv[0], N[0] + 1, endpoint=True)
        self.gridy = np.linspace(minv[1], maxv[1], N[1] + 1, endpoint=True)
        self.gridX = np.kron(self.gridx, np.ones(self.Ny + 1))
        self.gridY = np.kron(np.ones(self.Nx + 1), self.gridy)
        self.gridX = np.reshape(self.gridX, (self.Nx + 1, self.Ny + 1))
        self.gridY = np.reshape(self.gridY, (self.Nx + 1, self.Ny + 1))
        self.grid  = np.kron(self.gridy, self.gridx)
        self.grid  = np.reshape(self.grid, (self.Nx + 1, self.Ny + 1))
        self.grid  = np.flipud(self.grid)
        self.t = np.linspace(minv[2], maxv[2], N[2]+1, endpoint=True)

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print('plotting grid')

        plt.plot(self.gridX, self.gridY, 'b')
        plt.plot(self.gridY, self.gridX, 'b')
        plt.show()
