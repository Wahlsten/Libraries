def BackwardTimeCentralSpace(Nx, Nt, x, t, a, k, h, f, g):
    import numpy as np
    # A = [[b 0 0 ... 0 0 0],
    #     [c d 0 ... 0 0 0],
    #     [0 c d ... 0 0 0],
    #     ...
    #     [0 0 0 ... 0 c d],
    #
    # I = [e, 0 ... 0 h[1] 0 ... 0 h[2] 0 ... 0 h[-1]]

    b1 = -np.eye(Nx+1)


    b1[0, 0] = 0
    b1[0, 1] = 0
    b1[-1, -1] = 0
    b1[-1, -2] = 0

    b = np.diag(np.ones(Nt), -1)
    B = np.kron(b, b1)

    c = (1 + 2*a*k/(h*h)) * np.eye(Nx+1) - np.diag(np.ones(Nx), 1) * a * k / ( h * h) - np.diag(np.ones(Nx), -1) * a * k / ( h * h)

    c[0, 0] = 1
    c[0, 1] = 0
    c[-1, -1] = 1
    c[-1, -2] = 0

    c2 = np.diag(np.ones(Nt+1))
    c2[0,0] = 0
    C = np.kron(c2, c)

    d = np.zeros((Nt + 1, Nt + 1))

    d[0, 0] = 1
    D = np.kron(d, np.eye(Nx + 1))

    A = B + C + D

    e = np.zeros(Nt+1)
    e[0] = 1
    E = np.kron(e, f(x))

    h0 = np.zeros(Nx+1)
    h0[0] = 1
    H0 = np.kron(g(0, t), h0)
    H0[0] = 0

    h1 = np.zeros(Nx+1)
    h1[-1] = 1
    H1 = np.kron(g(1, t), h1)
    H1[Nx] = 0

    I = E + H0 + H1
    return A, I
