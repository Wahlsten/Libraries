def BackwardTimeBackwardSpace(Nx, Nt, x, t, a, k, h, f, g):
    import numpy as np
    # A = [[d 0 0 ... 0 0 0],
    #     [c b 0 ... 0 0 0],
    #     [0 c b ... 0 0 0],
    #     ...
    #     [0 0 0 ... 0 c b],
    #
    # I = [e, 0 ... 0 h[1] 0 ... 0 h[2] 0 ... 0 h[-1]]


    b = np.eye(Nx + 1)
    b[0, 0] = 0
    B = np.kron(b, np.eye(Nt + 1) - np.diag(np.ones(Nt), -1) * a * k / (h + a*k))

    c = -np.diag(np.ones(Nt+1)) * h / (h + a*k)
    c[0,0] = 0
    C = np.kron(np.diag(np.ones(Nx), -1), c)

    d = np.zeros((Nx + 1, Nx + 1))
    d[0, 0] = 1
    D = np.kron(d, np.eye(Nt + 1))

    A = B + C + D

    e = np.zeros(Nt+1)
    e[0] = 1
    E = np.kron(e, f(x))

    h = np.zeros(Nx+1)
    h[0] = 1
    H = np.kron(g(t), h)
    H[Nx] = 0

    I = E + H

    return A, I
