class gridClass:
    'Common base class for all employees'
    def __init__(self, minv, maxv, N):
        import numpy as np
        self.minx = minv[0]
        self.maxx = maxv[0]
        self.miny = minv[1]
        self.maxy = maxv[1]
        self.Lx   = N[0]
        self.Ly   = N[1]
        self.Nt   = N[2]

        self.generateMesh()

    def generateMesh(self):
        import numpy as np

        x = np.linspace(self.minx, self.maxx, self.Lx, endpoint=True)
        y = np.linspace(self.miny, self.maxy, self.Ly, endpoint=True)

        gridY, gridX = np.meshgrid(y, x)

        gridX = gridX.flatten()
        gridY = gridY.flatten()

        self.p = np.matrix([gridY, gridX])

        self.t = np.zeros((2*(self.Lx-1)*(self.Ly-1),4))
        k = 0
        for r in range(0,self.Ly-1):
            for c in range(0, self.Lx-1):

                self.t[k,   :] = [c + r*self.Lx,   c+1 + r*self.Lx,   (r+1)*self.Lx + c, 1]
                self.t[k+1, :] = [c+1 + r*self.Lx, (r+1)*self.Lx + c, (r+1)*self.Lx + c + 1, 1]
                k = k + 2

        self.t = np.transpose(self.t)
        A = range(0, self.Lx)
        A.extend(range(2*(self.Lx - 1)+1, self.Lx*self.Ly, self.Ly))
        A.extend(range(self.Lx * self.Ly-2, self.Lx * (self.Ly-1)-1, -1))
        A.extend(range(self.Lx*(self.Ly-2), 0, -self.Ly))

        B = A[1:]
        B.append(A[0])

        C = np.squeeze(np.asarray(np.zeros((1,len(B)))))
        D = np.squeeze(np.asarray(np.ones((1,len(B)))))
        E = range(0,len(B))
        F = np.squeeze(np.asarray(np.ones((1,len(B)))))
        G = np.squeeze(np.asarray(np.zeros((1,len(B)))))

        self.e = np.matrix([A,B,C,D,E,F,G])

        print self.t

    def plot_grid(self):
        import matplotlib.pyplot as plt
        import numpy as np
        print "plotting grid"

        numTri = len(np.squeeze(np.asarray(self.t[0, :])))
        X = np.zeros((4, numTri))
        Y = np.zeros((4, numTri))

        for k in range(0, numTri):

            X[0,k] = self.p[0, self.t[0, k]]
            X[1,k] = self.p[0, self.t[1, k]]
            X[2,k] = self.p[0, self.t[2, k]]
            X[3,k] = X[0,k]
            Y[0,k] = self.p[1, self.t[0, k]]
            Y[1,k] = self.p[1, self.t[1, k]]
            Y[2,k] = self.p[1, self.t[2, k]]
            Y[3,k] = Y[0,k]
            plt.plot(X[:,k], Y[:,k])

        xlen = (self.maxx - self.minx) * 0.1
        ylen = (self.maxy - self.miny) * 0.1
        plt.axis([self.minx - xlen, self.maxx + xlen, self.miny - ylen, self.maxy + ylen])
        plt.show()

class problemClass:
    def __init__(self, grid):
        import numpy as np

        def a_func(x): return (3 + 0*x)
        def b_func(x): return (1 + 0*x)

        self.a = a_func
        self.b = b_func

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
        print 'Initialize'

    def Scheme(self, grid, problem):
        print 'Create scheme'
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        import numpy as np

        N = len(np.squeeze(np.asarray(grid.p[0, :])))

        [S, M, b] = self.computeStiffnessMatrix(grid)

        bnodes2nodes = np.unique([grid.e[0,:], grid.e[1,:]])

        dofs = [x for x in range(1, N) if x not in bnodes2nodes]

        u_h = np.zeros((N,1))

        F = S[[[val] for val in dofs], [[dofs]]]
        F = F[0,:,:]
        B = np.squeeze(np.asarray(b[dofs]))

        u_h[dofs,0] = np.linalg.solve(F, B)

        print B

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_trisurf(np.squeeze(np.asarray(grid.p[0,:])), np.squeeze(np.asarray(grid.p[1,:])), np.squeeze(np.asarray(u_h.T)), linewidth=0.2, antialiased=True)
        plt.show()

    def boundaryConditionImposition(self):
        print 'Impose boundary conditions'

    def computeStiffnessMatrix(self, grid):
        print 'Compute stiffness matrix'
        import numpy as np

        def f(x): return x[0,:]*0 + x[1,:]*0 + 1

        M = len(np.squeeze(np.asarray(grid.t[0, :])))
        N = len(np.squeeze(np.asarray(grid.p[0,:])))

        St = np.zeros((N, N))
        Ma = np.zeros((N, N))
        b  = np.zeros((N, 1))

        [q, w] = self.quadrature(3, 'Triangle')

        [basis_q, grad_basis_q] = self.basisLinear(q)

        local_M_mat = np.zeros((3, 3))
        local_S_mat = np.zeros((3, 3))
        local_load = np.zeros((3, 1))

        u = np.squeeze(np.asarray(np.zeros((1, 5))))

        for k in range(0, M-1):

            nodes = grid.p[:, np.asarray(grid.t[:3, k], dtype = int)]

            A = np.squeeze(np.asarray(nodes[:, 1] - nodes[:, 0]))
            B = np.squeeze(np.asarray(nodes[:, 2] - nodes[:, 0]))

            C = np.matrix([A, B])
            cv = nodes[:, 0]
            Cdet = C[0, 0] * C[1, 1] - C[0, 1] * C[1, 0]

            CinvT = 1.0/Cdet * np.array([[C[1, 1], -C[1, 0]], [-C[0, 1], C[0, 0]]])

            q_loc = np.dot(C, q) + np.kron(cv, np.ones((1, len(np.squeeze(np.asarray(q[0,:]))))))

            f_at_q = np.squeeze(np.asarray(f(q_loc)))

            for i in range(0, 2):

                local_load[i] = np.dot((f_at_q * basis_q[i, :]), w) * Cdet

                for j in range(0, 2):

                    local_M_mat[i, j] = np.dot(basis_q[i,:] * basis_q[j,:], w) * Cdet
                    local_S_mat[i, j] = np.dot(np.dot(CinvT, grad_basis_q[:,0,i]), \
                    np.dot(CinvT, grad_basis_q[:,0,j])) * Cdet


            tp = np.squeeze(np.asarray(grid.t[0:3, k], dtype = int))
            b[tp] = b[tp] + local_load
            St[[[val] for val in tp], [[tp]]] = St[[[val] for val in tp], [[tp]]] + local_S_mat
            Ma[[[val] for val in tp], [[tp]]] = Ma[[[val] for val in tp], [[tp]]] + local_M_mat

        print b
        return St, Ma, b

    def basisLinear(self, x):
        print 'Compute linear basis'
        import numpy as np

        M = len(np.squeeze(np.asarray(x[0,:])))

        value = np.zeros((3, M))

        value[0, :] = np.ones((1,M)) - x[0, :] - x[1, :]
        value[1, :] = x[0, :]
        value[2, :] = x[1, :]

        d_value = np.zeros((2, M, 3))
        v = np.ones((1, M))

        d_value[:,:,0] = [-v, -v]
        d_value[:,:,1] = [v, np.zeros((1,M))]
        d_value[:,:,2] = [np.zeros((1,M)), v]

        return value, d_value

    def quadrature(self, order, elementType):
        import numpy as np

        if elementType is 'Triangle':
            return self.quadTriangle(order)
        elif elementType is 'Quadrilatural':
            return self.quadQuadlirateral(order)

    def quadTriangle(self, order):
        import numpy as np

        if order is 1:
            q = np.matrix([[1/3., 1/3.]])
            w = np.array([1/2.])
        elif order is 3:
            q = np.transpose(np.matrix([[0.5, 0], [0, 0.5], [0.5, 0.5]]))
            w = np.array([1/6., 1/6., 1/6.])
        elif order is 4:
            q = np.matrix([[1/3., 1/3.], [0.6, 0.2], [0.2, 0.6], [0.2, 0.2]])
            w = np.array([-27/96., 25/96., 25/96., 25/96])

        return q, w

    def quadQuadlirateral(self):
        import numpy as np

        if order is 1:
            q = np.matrix([[0, 0]])
            w = np.array([2])
        elif order is 4:
            q = np.matrix([[-np.sqrt(1/3.), -np.sqrt(1/3.)], [np.sqrt(1/3.), -np.sqrt(1/3.)], [-np.sqrt(1/3.), np.sqrt(1/3.)], [np.sqrt(1/3.), np.sqrt(1/3.)]])
            w = np.array([1, 1, 1, 1])
        elif order is 9:
            q = np.matrix([[-np.sqrt(3/5.), -np.sqrt(3/5.)], [0, -np.sqrt(3/5.)], [np.sqrt(3/5.), -np.sqrt(3/5.)], [-np.sqrt(3/5.), 0], [0, 0], \
            [np.sqrt(3/5.), 0], [-np.sqrt(3/5.), np.sqrt(3/5.)], [0, np.sqrt(3/5.)], [np.sqrt(3/5.), np.sqrt(3/5.)]])
            w = np.array([25/81., 40/81., 25/81., 40/81., 64/81., 40/81., 25/81., 40/81., 25/81.])

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

g = gridClass([0, 0, 0], [1, 1, 1], [10, 10, 4000])
g.plot_grid()

p = problemClass(g)
p.createData(g)

s = SchemeClass(p)
s.Scheme(g, p)
s.Error(g, p)
