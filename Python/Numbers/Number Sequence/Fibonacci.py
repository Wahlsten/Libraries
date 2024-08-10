import numpy as np
def Fibonacci(N):

    Numbers = np.zeros(N)
    Numbers[0] = 1
    Numbers[1] = 1

    for n in range(2, N):

        Numbers[n] = Numbers[n-1] + Numbers[n-2]

    return Numbers
