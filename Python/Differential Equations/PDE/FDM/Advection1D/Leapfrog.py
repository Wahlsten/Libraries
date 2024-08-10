def Leapfrog(Nx, Nt, x, t, a, k, h, f, g):
    import numpy as np
    # A = [[b 0 0 ... 0 0 0],
    #     [c d 0 ... 0 0 0],
    #     [0 c d ... 0 0 0],
    #     ...
    #     [0 0 0 ... 0 c d],
    #
    # I = [e, 0 ... 0 h[1] 0 ... 0 h[2] 0 ... 0 h[-1]]

    b1 = np.diag(np.ones(Nt), 1) * a * k / h - np.diag(np.ones(Nt), -1) * a * k / h

    if a > 0:
        b1[0, 1] = 0
        b1[-1, -2] = -2*a*k/h
        b1[-1, -1] =  2*a*k/h
    elif a < 0:
        b1[-1, -2] = 0
        b1[0, 1] =  2*a*k/h
        b1[0, 0] = -2*a*k/h

    b = np.eye(Nx + 1)
    b[0, 0] = 0
    b[-1, -1] = 0
    B = np.kron(b, b1)

    b3 = np.eye(Nt+1) + np.diag(np.ones(Nt), 1) * a * k / (2*h) - np.diag(np.ones(Nt), -1) * a * k / (2*h)

    if a < 0:
        b3[0,0] = 1 - a*k/h
        b3[-1,-2] = 0
    elif a > 0:
        b3[-1,-1] = 1 + a*k/h
        b3[0,1] = 0
    b33 = np.zeros((Nx + 1, Nx + 1))
    b33[-1,-1] = 1

    B3 = np.kron(b33, b3)

    b2 = np.eye(Nt + 1)

    if a > 0:
        b2[0, 0] = 0

    elif a < 0:
        b2[-1, -1] = 0

    b22 = np.ones(Nx)
    b22[0] = 0
    b22 = np.diag(b22, 1)
    B2 = np.kron(b22, b2)

    c = -np.diag(np.ones(Nt+1))

    if a > 0:
        c[0, 0] = 0
    elif a < 0:
        c[-1,-1] = 0

    C = np.kron(np.diag(np.ones(Nx), -1), c)

    d = np.zeros((Nx + 1, Nx + 1))

    d[0, 0] = 1
    D = np.kron(d, np.eye(Nt + 1))

    b4 = np.eye(Nt+1)
    b4[0,0] = 0
    b4[-1,-1] = 0
    B4 = np.zeros((Nx+1, Nx+1))
    if a < 0:
        B4[-1,-1] = 1
    elif a > 0:
        B4[0, 0] = 1
    B4 = np.kron(b4, B4)

    A = B + B2 + B3 + B4 + C + D

    e = np.zeros(Nt+1)
    e[0] = 1
    E = np.kron(e, f(x))

    h = np.zeros(Nx+1)
    if a > 0:
        h[0] = 1
    elif a < 0:
        h[-1] = 1
    H = np.kron(g(t), h)

    if a > 0:
        H[0] = 0
    elif a < 0:
        H[Nx] = 0

    I = E + H

    return A, I
