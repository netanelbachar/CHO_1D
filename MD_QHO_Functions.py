import numpy as np
import matplotlib.pyplot as plt
import math
import time


def initial_position(mass, w, hbar, beads):
    mu, sigma = 0, math.sqrt(hbar / (mass * w))
    x = np.random.normal(mu, sigma, beads)
    return x


def initial_velocity(mass, beta, beads):
    mu, sigma = 0, math.sqrt(1 / (beta * mass))
    vx = np.random.normal(mu, sigma, beads)
    return vx


def oscillator_force(mass, w, beads, x):
    f1 = - mass * w**2 * x / beads
    return f1


def ring_springs_force(mass, wp, x):
    c2 = - mass * wp**2
    f2 = np.zeros(len(x))
    for j in range(len(x)):
        if j == (len(x) - 1):
            f2[j] = c2 * (2 * x[j] - x[0] - x[j - 1])
        elif j == 0:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[-1])
        else:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[j - 1])
    return f2


def oscillator_potential(mass, w, x):
    potential = (0.5 * mass * w**2 * x**2).sum()
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


def xsi(beads):  # This is the noise term
    z = np.random.normal(0, 1, beads)
    return z


def langevin(mass, beta, v, dt,beads):
    # is xsi constant for entire time propagation?
    c1 = math.exp(-1 * gamma(dt) * dt / 2)
    c2 = math.sqrt((1 / beta / mass) * (1 - c1**2))
    vel = c1 * v + c2 * xsi(beads)
    return vel


def kinetic_estimator(beta, beads, mass, wp, x):
    # fac = beads / 5
    k2 = 0
    c2 = 0.5 * mass * wp**2
    for j in range(len(x)):
        if j == (len(x) - 1):
            k2 += c2 * (x[j] - x[0])**2
        else:
            k2 += c2 * (x[j] - x[j+1])**2
    k_estimator = (beads / (2 * beta)) - k2
    return k_estimator  # / fac


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
    e_change = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    pot_est = np.zeros(g_steps)
    kin_est = np.zeros(g_steps)
    # Histogram
    pos = np.zeros(g_steps)
    vel = np.zeros(g_steps)

    for step in range(0, g_steps):
        times[step] = t
        steps[step] = int(s)
        kin[step] = kinetic_1d(mass, vx)
        potential[step] = oscillator_potential(mass, w, x) / beads + springs_potential(mass, wp, x)
        e_tot[step] = kin[step] + potential[step]
        e_change[step] = abs(e_tot[step] - e_tot[0]) * 100 / e_tot[0]
        temp_exp[step] = 2.0 * kin[step] / (kboltz * beads)  # Divided by beads from the Equipartition Function
        pot_est[step] = oscillator_potential(mass, w, x) / beads
        kin_est[step] = kinetic_estimator(beta, beads, mass, wp, x)
        t += dt
        s += 1
        # Histogram
        pos[step] = x[0]
        vel[step] = vx[0]
        
        # Time propagation
        vx = langevin(mass, beta, vx, dt, beads)
        vx = vx + 0.5 * dt * (force / mass)
        x = x + vx * dt
        force = oscillator_force(mass, w, beads, x) + ring_springs_force(mass, wp, x)
        vx = vx + 0.5 * dt * (force / mass)
        vx = langevin(mass, beta, vx, dt, beads)

    return steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est



