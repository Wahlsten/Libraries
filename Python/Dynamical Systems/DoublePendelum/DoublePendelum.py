import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# Set up the figure 
# the axis
# the plot element
fig = plt.figure()
ax = plt.axes(xlim=(-400, 400), ylim=(-400, 400))
line2, = ax.plot([], [], 'k-', lw=2)
line3, = ax.plot([], [], 'b-', lw=2)
line4, = ax.plot([], [], 'b-', lw=2)
line5, = ax.plot([], [], 'b-', lw=2)
line6, = ax.plot([], [], 'b-', lw=2)

# Initialization function:
# Plot the background of each frame
def init():
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    line6.set_data([], [])
    return line2, line3, line4, line5, line6

def PendulumPosition(dt):

    global theta_1
    global vel_1
    global acc_1
    global theta_2
    global vel_2
    global acc_2
    global length_1
    global length_2
    global mass_1
    global mass_2
    global g
    global prev_x
    global prev_y

    if dt == 0:
        length_1 = 200.
        length_2 = 200.
        mass_1 = 40.
        mass_2 = 40.
        g = 1.
        theta_1 = np.pi/1.65
        theta_2 = np.pi/1.65
        vel_1 = 0.
        acc_1 = 0.
        acc_2 = 0.
        vel_2 = 0.
        prev_x = np.zeros(2000)
        prev_y = np.zeros(2000)

    else:

        num1 = - g * (2. * mass_1  + mass_2) * np.sin(theta_1)
        num2 = - mass_2 * g * np.sin(theta_1 - 2. * theta_2)
        num3 = - 2. * np.sin(theta_1 - theta_2) * mass_2
        num4 = (vel_2 * vel_2 * length_2 + vel_1 * vel_1 * length_1 * np.cos(theta_1 - theta_2))
        den  = length_1 * (2. * mass_1 + mass_2 - mass_2 * np.cos(2 * theta_1 - 2. * theta_2))

        num1_2 = 2. * np.sin(theta_1 - theta_2)
        num2_2 = vel_1 * vel_1 * length_1 * (mass_1 + mass_2)
        num3_2 = g * (mass_1 + mass_2) * np.cos(theta_1)
        num4_2 = vel_2 * vel_2 * length_2 * mass_2 * np.cos(theta_1 - theta_2)
        den_2  = length_2 * (2. * mass_1 + mass_2 - mass_2 * np.cos(2. * theta_1 - 2. * theta_2))

        acc_1 = (num1 + num2 + num3 * num4) / den
        acc_2 = (num1_2 * (num2_2 + num3_2 + num4_2)) / den_2
        vel_1 = vel_1 + acc_1
        vel_2 = vel_2 + acc_2
        theta_1 = theta_1 + vel_1
        theta_2 = theta_2 + vel_2

    return theta_1, theta_2

def PendelumDraw(i, theta_1, theta_2):

    global length_1
    global length_2
    global mass_1
    global mass_2
    global prev_x
    global prev_y

    draw_mass_1 = mass_1/3.
    draw_mass_2 = mass_2/3.

    x_pos_1 = np.sin(theta_1) * length_1
    x_pos_2 = np.sin(theta_2) * length_2

    y_pos_1 = -np.cos(theta_1) * length_1
    y_pos_2 = -np.cos(theta_2) * length_2

    prev_x[i] = x_pos_1 + x_pos_2
    prev_y[i] = y_pos_1 + y_pos_2
    x_p = list(reversed(prev_x[0:i+1]))
    y_p = list(reversed(prev_y[0:i+1]))
    x = [0, x_pos_1, x_pos_1 + x_pos_2]
    y = [0, y_pos_1, y_pos_1 + y_pos_2]

    x_l1 = [0, x_pos_1 - np.sin(theta_1) * draw_mass_1]
    y_l1 = [0, y_pos_1 + np.cos(theta_1) * draw_mass_1]

    x_l2 = [x_pos_1 + np.sin(theta_2) * draw_mass_1,  x_pos_1 + x_pos_2 - np.sin(theta_2) * draw_mass_2]
    y_l2 = [y_pos_1 - np.cos(theta_2) * draw_mass_1,  y_pos_1 + y_pos_2 + np.cos(theta_2) * draw_mass_2]

    t_size = 100
    t = np.linspace(0, 2*np.pi, t_size)
    x_c1 = draw_mass_1 * np.sin(t) + np.ones(t_size) * x_pos_1
    y_c1 = draw_mass_1 * np.cos(t) + np.ones(t_size) * y_pos_1

    x_c2 = draw_mass_2 * np.sin(t) + np.ones(t_size) * (x_pos_1 + x_pos_2)
    y_c2 = draw_mass_2 * np.cos(t) + np.ones(t_size) * (y_pos_1 + y_pos_2)

    return x_l1, y_l1, x_l2, y_l2, x_p, y_p, x_c1, y_c1, x_c2, y_c2

# Animation function
# This is called sequentially
def animate(i):

    theta_1, theta_2 = PendulumPosition(i)
    x_l1, y_l1, x_l2, y_l2, x_p, y_p, x_c1, y_c1, x_c2, y_c2 = PendelumDraw(i, theta_1, theta_2)
    #x = np.linspace(0, 2, 1000)
    #y = np.sin(2 * np.pi * (x - 0.01 * i))

    line2.set_data(x_p, y_p)
    line3.set_data(x_c1, y_c1)
    line4.set_data(x_c2, y_c2)
    line5.set_data(x_l1, y_l1)
    line6.set_data(x_l2, y_l2)
    return line2, line3, line4, line5, line6

# Call the animator. 
# blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2000, interval=20, blit=True)


plt.show()
