import numpy as np
import matplotlib.pyplot as plt
from RungeKutta4 import RungeKutta4

# parameters
beta  = 1.0          # expected amount of people an infected person infects per day 
D     = 4.0          # number of days an infected person has and can spread the disease
gamma = 1.0 / D      # the proportion of infected recovering per day
R_0   = beta / gamma # the total number of people an infected person infects
N     = 1000         # total population

def ft(t, y):
    S = y[0]
    I = y[1]
    R = y[2]

    S_eq = -beta * I * S / N
    I_eq = beta * I * S / N - gamma * I
    R_eq = gamma * I
    
    return np.array([S_eq, I_eq, R_eq])

# Initial conditions
S_0 = 999 # number of people susceptible on day 0
I_0 = 1   # number of people infected on day 0
R_0 = 0   # number of people recovered on day 0

num_gridpoints = 1000 # number of grid points
T              = 40 # num days

[odeT, odeY] = RungeKutta4([0, T], [S_0, I_0, R_0], ft, num_gridpoints)

plt.plot(odeT, odeY[0,:])
plt.plot(odeT, odeY[1,:])
plt.plot(odeT, odeY[2,:])
plt.legend(['Susceptible', 'Infected', 'Recovered'])
plt.show()

