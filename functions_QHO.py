import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import time


def potential_classic(mass, w, x):
    potential = 0.5 * mass * w**2 * x**2
    return potential


def kinetic_classic(mass, vx):
    kinetic_x = 0
    kinetic_x = 0.5 * mass * vx ** 2
    return kinetic_x


def force_classic(mass, w, x):
    f = - mass * w**2 * x
    return f

def verlet_classic(g_steps, mass, dt, beta, w, hbar, kboltz, beads):
    t = 0
    x = initial_position(mass, w, hbar, beads)
    vx = velocity(mass, beta, beads)
    force_x = force_classic(mass, w, x)
    times = np.zeros(g_steps)
    pos = np.zeros(g_steps)
    vx_n = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    pot = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)

    for step in range(0, g_steps):
        times[step] = t
        t += dt
        kin[step] = kinetic_classic(mass, vx)
        pot[step] = potential_classic(mass, w, x)
        e_tot[step] = kin[step] + pot[step]
        pos[step] = x
        vx_n[step] = vx
        temp_exp[step] = kboltz * e_tot[step]

        # Verlet
        vx = vx + 0.5 * dt * force_x / mass
        x = x + dt * vx
        force_x = force_classic(mass, w, x)
        vx = vx + 0.5 * dt * force_x / mass

    return times, pos, vx_n, kin, pot, e_tot, temp_exp


def langevin_classic(g_steps, mass, dt, beta, w, hbar, kboltz, beads):
    t = 0
    x = initial_position(mass, w, hbar, beads)
    vx = velocity(mass, beta, beads)
    force_x = force_classic(mass, w, x)
    times = np.zeros(g_steps)
    pos = np.zeros(g_steps)
    vx_n = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    pot = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)

    for step in range(0, g_steps):
        times[step] = t
        t += dt
        kin[step] = kinetic_classic(mass, vx)
        pot[step] = potential_classic(mass, w, x)
        e_tot[step] = kin[step] + pot[step]
        pos[step] = x
        vx_n[step] = vx
        temp_exp[step] = kboltz * e_tot[step]

        # Verlet
        vx = langevin(mass, beta, vx, dt)
        vx = vx + 0.5 * dt * force_x / mass
        x = x + dt * vx
        force_x = force_classic(mass, w, x)
        vx = vx + 0.5 * dt * force_x / mass
        vx = langevin(mass, beta, vx, dt)

    return times, pos, vx_n, kin, pot, e_tot, temp_exp


def initial_position(mass, w, hbar, beads):
    """

    :return: Position Array with gaussian distribution
    length as the number of beads
    """
    mu = 0
    sigma = math.sqrt(hbar / (mass * w))
    x = np.random.normal(mu, sigma, beads)
    # x = np.zeros(beads)
    # x[0] = 0
    return x


def velocity(mass, beta, beads):
    """

    :return: Velocity Array of zeros.
    length as the number of beads
    """
    mu = 0
    sigma = math.sqrt(1 / (beta * mass))
    vx = np.random.normal(mu, sigma, beads)
    # vx = np.zeros(beads)
    # vx[0] = 0

    return vx


def force(mass, w, wp, beads, x):
    """

    :param mass: mass of partilce
    :param w: quantum oscillator freq
    :param wp: spring between bead freq
    :param beads: number of beads
    :param x: array of position
    :return:
    """
    f2 = np.zeros(len(x))
    f1 = - mass * w ** 2 * x / beads
    for j in range(len(x)):
        if j == (len(x) - 1):
            f2[j] = - mass * wp**2 * (2 * x[j] - x[0] - x[j - 1])
        elif j == 0:
            f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[-1])
        else:
            f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[j - 1])
    force_beads = f2 + f1
    return force_beads


def potential_1d(mass, w, wp, beads, x):
    """

    :param x: Array of bead positions
    :return: Sum of all potential energy of the beads
    """
    potential = 0
    for j in range(len(x)):
        if j == (len(x) - 1):
            potential += (0.5 * mass * wp ** 2 * (x[0] - x[j])**2) + ((mass * w**2) / (2 * beads)) * x[j] ** 2
        else:
            potential += (0.5 * mass * wp ** 2 * (x[j + 1] - x[j])**2) + ((mass * w**2) / (2 * beads)) * x[j] ** 2
    return potential


def kinetic_1d(mass, v):
    """

    :param v: Velocity array of the beads
    :return: Total Kinetic energy of the beads
    """
    kinetic = (0.5 * mass * v**2).sum()
    # kinetic = 0
    # for j in range(len(v)):
    #     kinetic += 0.5 * mass * v[j] ** 2
    return kinetic


# Langevin Dynamics
def l_gamma(dt):  # This is the friction term
    gamma = 0.01 / dt
    return gamma


def l_xsi():  # This is the noise term
    xsi = np.random.normal(0, 1, 1)
    return xsi


def langevin(mass, beta, v, time_step):
    c1 = math.exp(-1 * l_gamma(time_step) * time_step / 2)
    c2 = math.sqrt((mass / beta) * (1 - c1**2))
    vel = c1 * v + c2 * l_xsi()
    return vel




def verlet_propagator(g_steps, dt, mass, w, beta, hbar, kboltz, wp, beads):
    """

    :param g_steps: number of time steps
    :param dt: time step
    :param mass: mass of bead
    :param w: oscillator freq
    :param beta: inverse temperature
    :param hbar: plack constant
    :param kboltz: boltzman constant
    :param wp: between beads freq
    :param beads: number of beads
    :return: steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp
    """
    t = 0
    s = 0
    x = initial_position(mass, w, hbar, beads)
    vx = velocity(mass, beta, beads)
    force_x = force(mass, w, wp, beads, x)

    times = np.zeros(g_steps)
    steps = np.zeros(g_steps)
    pos = np.zeros(g_steps)
    vx_n = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    pot = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    percent_change = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    for step in range(0, g_steps):
        steps[step] = s
        times[step] = t
        kin[step] = kinetic_1d(mass, vx)
        pot[step] = potential_1d(mass, w, wp, beads, x)
        e_tot[step] = kin[step] + pot[step]
        pos[step] = x[0]
        vx_n[step] = vx[0]
        percent_change[step] = abs(e_tot[step] - e_tot[0]) * 100 / e_tot[0]
        temp_exp[step] = kboltz * e_tot[step] / beads  # Divided by beads from the Equipartition Function
        s += 1
        t += dt
        # Time propagation
        vx = vx + 0.5 * dt * force_x / mass
        x = x + dt * vx
        force_x = force(mass, w, wp, beads, x)
        vx = vx + 0.5 * dt * force_x / mass

    return steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp



def langevin_thermostat(g_steps, dt, mass, w, beta, hbar, kboltz, wp, beads):
    """

    :param g_steps: number of time steps
    :param dt: time step
    :param mass: mass of bead
    :param w: oscillator freq
    :param beta: inverse temperature
    :param hbar: plack constant
    :param kboltz: boltzman constant
    :param wp: between beads freq
    :param beads: number of beads
    :return: time array, pos array only of bead number 1, vx_n velocity of bead number 1
    , total energy (kinetic and potential) of entire ring, temperature of system using equipartition function
    , potential estimator for a given time t.
    """
    t = 0
    s = 0
    x = initial_position(mass, w, hbar, beads)
    vx = velocity(mass, beta, beads)
    force_x = force(mass, w, wp, beads, x)

    times = np.zeros(g_steps)
    steps = np.zeros(g_steps)
    pos = np.zeros(g_steps)
    vx_n = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    pot = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    percent_change = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    pot_est = np.zeros(g_steps)
    kin_est = np.zeros(g_steps)

    for step in range(0, g_steps):
        steps[step] = s
        times[step] = t
        kin[step] = kinetic_1d(mass, vx)
        pot[step] = potential_1d(mass, w, wp, beads, x)
        e_tot[step] = kin[step] + pot[step]
        pos[step] = x[0]
        vx_n[step] = vx[0]
        percent_change[step] = abs(e_tot[step] - e_tot[0]) * 100 / e_tot[0]
        temp_exp[step] = kboltz * e_tot[step] / beads  # Divided by beads from the Equipartition Function
        pot_est[step] = pot_estimator(mass, w, beads, x)
        kin_est[step] = kinetic_estimator(beta, beads, mass, wp, x)
        s += 1
        t += dt
        # Time propagation
        vx = langevin(mass, beta, vx, dt)
        vx = vx + 0.5 * dt * force_x / mass
        x = x + dt * vx
        force_x = force(mass, w, wp, beads, x)
        vx = vx + 0.5 * dt * force_x / mass
        vx = langevin(mass, beta, vx, dt)

    return steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp, pot_est, kin_est

 # times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp, pot_est

def pot_estimator(mass, w, beads, x):
    '''
    :param mass: mass of particle
    :param w: Oscilator freq
    :param beads: number of beads
    :param x: array of bead positions
    :return: Potential Estimator in for a given t.
    '''
    estimator = sum(0.5 * mass * w**2 * x**2) / beads
    return estimator


def estimator_potential_array(g_steps, dt, mass, w, beta, hbar, kboltz, beads_array):
    pot_est_array = np.zeros(len(beads_array))
    for j, p in enumerate(beads_array):  # j index and p is the number of beads
        wp = math.sqrt(p) / (beta * hbar)
        steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp, pot_est, kin_est = \
            langevin_thermostat(g_steps, dt, mass, w, beta, hbar, kboltz, wp, p)
        # Throw g_step*0.1 first variables
        pot_est_avg = sum(pot_est[int(g_steps*0.4):])/len(pot_est[int(g_steps*0.4):])
        pot_est_array[j] = pot_est_avg
    return pot_est_array


def kinetic_estimator(beta, beads, mass, wp, x):
    for j in range(len(x)):
        if j == (len(x) - 1):
            k_estimator = (beads / (2 * beta)) - 0.5 * mass * wp**2 * (x[0]-x[j])**2
        else:
            k_estimator = (beads / (2 * beta)) - 0.5 * mass * wp**2 * (x[j+1]-x[j])**2
    return k_estimator







## DONt touch

# def oscillator_force(mass, w, beads, x):
#     f1 = - mass * w**2 * x / beads
#     return f1
#
#
# def ring_springs_force(mass, wp, x):
#     f2 = np.zeros(len(x))
#     for j in range(len(x)):
#         if j == (len(x) - 1):
#             f2[j] = - mass * wp**2 * (2 * x[j] - x[0] - x[j - 1])
#         elif j == 0:
#             f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[-1])
#         else:
#             f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[j - 1])
#     return f2

