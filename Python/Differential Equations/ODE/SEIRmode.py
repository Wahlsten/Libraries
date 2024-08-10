import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters
alpha = 0.2          # Death rate
beta  = 1.0          # expected amount of people an infected person infects per day 
delta = 1.0 / 3.0    # length of incubation period
rho   = 1.0 / 9.0    # rate at which people die (1/days from infected until death)
D     = 4.0          # number of days an infected person has and can spread the disease
gamma = 1.0 / D      # the proportion of infected recovering per day
R_0   = beta / gamma # the total number of people an infected person infects
N     = 1000000      # total population

def ft(y, t, N, beta, gamma, delta, alpha, rho):
    S, E, I, R, D = y

    dS_dt = -beta * I * S / N
    dE_dt = beta * I * S / N - delta * E
    dI_dt = delta * E - (1 - alpha) * gamma * I - alpha * rho * I
    dR_dt = (1 - alpha) * gamma * I
    dD_dt = alpha * rho * I
    
    return np.array([dS_dt, dE_dt, dI_dt, dR_dt, dD_dt])

# Initial conditions
S_0 = N - 1 # number of people susceptible on day 0
I_0 = 1     # number of people infected on day 0
E_0 = 0     # number of exposed people on day 0
R_0 = 0     # number of people recovered on day 0
D_0 = 0     # number of dead people on day 0

num_gridpoints = 100 # number of grid points
T              = 100 # num days
t              = np.linspace(0, T, num_gridpoints)

y_0 = S_0, E_0, I_0, R_0, D_0

solution = odeint(ft, y_0, t, args=(N, beta, gamma, delta, alpha, rho))

S, E, I, R, D = solution.T

plt.plot(t, S)
plt.plot(t, E)
plt.plot(t, I)
plt.plot(t, R)
plt.plot(t, D)
plt.legend(['Susceptible', 'Exposed', 'Infected', 'Recovered', 'Deaths'])
plt.show()

