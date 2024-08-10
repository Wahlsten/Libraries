import numpy as np

def NelderMead(f, tol, N):

    Nmin = 1
    Nmax = N

    V = np.array([[0., 0.], [0.1, 0.1], [1., -0.1]])
    [mm, n] = np.shape(V)

    Y = np.zeros(n+1)
    P = np.zeros((Nmax, n))
    Q = np.zeros(Nmax)

    for k in range(0, n+1):

        Y[k] = f(V[k,0], V[k,1])

    mn = np.min(Y)
    mx = np.max(Y)
    lo = [i for i, j in enumerate(Y) if j == mn]
    hi = [i for i, j in enumerate(Y) if j == mx]

    lo = lo[0]
    hi = hi[0]
    li = hi + 0
    ho = lo + 0

    for j in range(0, n + 1):

        if ((j != lo) & (j != hi)) & (Y[j] <= Y[li]):
            li = j

        if ((j != hi) & (j != lo)) & (Y[j] >= Y[ho]):
            ho = j

    cnt = 0

    while ((Y[hi] > Y[lo] + tol) & (cnt < Nmax)) | (cnt < Nmin):

        S = np.zeros(n)

        for j in range(0, n + 1):

            S = S + V[j,:]

        M = (S - V[hi, :])/float(n)
        R = 2*M - V[hi, :]
        yR = f(R[0], R[1])

        if (yR < Y[ho]):

            if Y[li] < yR:

                V[hi, :] = R
                Y[hi] = yR

            else:

                E = 2*R - M
                yE = f(E[0], E[1])

                if yE < Y[li]:

                    V[hi, :] = E
                    Y[hi] = yE

                else:

                    V[hi, :] = R
                    Y[hi] = yR

        else:

            if (yR < Y[hi]):

                V[hi, :] = R
                Y[hi] = yR

            C = (V[hi, :] + M) / 2.
            yC = f(C[0], C[1])
            C2 = (M + R)/2.
            yC2 = f(C2[0], C2[1])

            if yC2 < yC:

                V[hi, :] = C
                Y[hi] = yC

            else:

                for j in range(0, n + 1):

                    if (j != lo):

                        V[j, :] = (V[j, :] + V[lo, :]) / 2.
                        Z = V[j, :]
                        Y[j] = f(Z[0], Z[1])

        mn = np.min(Y)
        mx = np.max(Y)

        lo = [i for i, j in enumerate(Y) if j == mn]
        hi = [i for i, j in enumerate(Y) if j == mx]

        lo = lo[0]
        hi = hi[0]

        li = hi + 0
        ho = lo + 0

        for j in range(0, n + 1):

            if ((j != lo) & (j != hi)) & (Y[j] <= Y[li]):
                li = j

            if ((j != hi) & (j != lo)) & (Y[j] >= Y[ho]):
                ho = j

        P[cnt, :] = V[lo, :]
        Q[cnt] = Y[lo]

        cnt = cnt + 1


    xmin = P[cnt-1,:]
    minimum  = Q[cnt-1]

    return minimum, xmin
