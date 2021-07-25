import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import time


def initial_position(mass, w, hbar, beads):
    mu = 0
    sigma = math.sqrt(hbar / (mass * w))
    x = np.random.normal(mu, sigma, beads)
    return x


def initial_velocity(mass, beta, beads):
    mu = 0
    sigma = math.sqrt(1 / (beta * mass))
    vx = np.random.normal(mu, sigma, beads)
    return vx


def oscillator_force(mass, w, beads, x):
    f1 = - mass * w**2 * x / beads
    return f1


def ring_springs_force(mass, wp, x):
    f2 = np.zeros(len(x))
    for j in range(len(x)):
        if j == (len(x) - 1):
            f2[j] = - mass * wp**2 * (2 * x[j] - x[0] - x[j - 1])
        elif j == 0:
            f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[-1])
        else:
            f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[j - 1])
    return f2

def trying(mass, wp, x):
    f2 = np.zeros(len(x))
    for j in range(0, len(x)):
        f2[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[j - 1])



def oscillator_potential(mass, w, x):
    potential = (((mass * w**2) / 2) * x**2).sum()
    return potential


def springs_potential(mass, wp, x):
    potential = 0
    for j in range(len(x)):
        if j == (len(x) - 1):
            potential += 0.5 * mass * wp**2 * (x[0] - x[j])**2
        else:
            potential += 0.5 * mass * wp**2 * (x[j + 1] - x[j])**2
    return potential


def kinetic_1d(mass, v):
    kinetic = (0.5 * mass * v**2).sum()
    return kinetic


def gamma(dt):  # This is the friction term
    g = 1 / (100 * dt)
    return g


def xsi():  # This is the noise term
    z = np.random.normal(0, 1, 1)
    return z


def langevin(mass, beta, v, dt):
    c1 = math.exp(-1 * gamma(dt) * dt / 2)
    c2 = math.sqrt((mass / beta) * (1 - c1**2))
    vel = c1 * v + c2 * xsi()
    return vel


def potential_estimator(mass, w, beads, x):
    estimator = (0.5 * mass * w**2 * x**2).sum()
    estimator = estimator / beads
    return estimator


def kinetic_estimator(beta, beads, mass, wp, x):
    k_estimator = 0
    for j in range(len(x)):
        if j == (len(x) - 1):
            k_estimator += (beads / (2 * beta)) - 0.5 * mass * wp**2 * (x[0]-x[j])**2
        else:
            k_estimator += (beads / (2 * beta)) - 0.5 * mass * wp**2 * (x[j+1]-x[j])**2
    return k_estimator


def langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads):
    wp = math.sqrt(beads) / (beta * hbar)
    x = initial_position(mass, w, hbar, beads)
    vx = initial_velocity(mass, beta, beads)
    force = oscillator_force(mass, w, beads, x) + ring_springs_force(mass, wp, x)

    t = 0
    s = 0
    times = np.zeros(g_steps)
    steps = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    potential = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    pot_est = np.zeros(g_steps)
    kin_est = np.zeros(g_steps)

    for step in range(0, g_steps):
        times[step] = t
        steps[step] = int(s)
        kin[step] = kinetic_1d(mass, vx)
        potential[step] = oscillator_potential(mass, w, x) / beads + springs_potential(mass, wp, x)
        e_tot[step] = kin[step] + potential[step]
        temp_exp[step] = e_tot[step] / (kboltz * beads)  # Divided by beads from the Equipartition Function
        pot_est[step] = potential_estimator(mass, w, beads, x)
        kin_est[step] = kinetic_estimator(beta, beads, mass, wp, x)
        t += dt
        s += 1

        # Time propagation
        vx = langevin(mass, beta, vx, dt)
        vx = vx + 0.5 * dt * (force / mass)
        x = x + vx * dt
        force = oscillator_force(mass, w, beads, x) + ring_springs_force(mass, wp, x)
        vx = vx + 0.5 * dt * (force / mass)
        vx = langevin(mass, beta, vx, dt)

    return steps, times, kin, e_tot, temp_exp, pot_est, kin_est








