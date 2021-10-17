import numpy as np
import matplotlib.pyplot as plt
import math
import time

mass, w, hbar, kboltz = 1, 1, 1, 1
beta = 0.3

N_particles = 2
# beads = 32

# Time
T_max = 40000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)

# ErrorBars
cutoff = int(g_steps * 0.4)
num_blocks = 120
block_size = int((g_steps-cutoff)/num_blocks)


seed = 347885
np.random.seed(seed)
# wp = math.sqrt(beads) / (beta * hbar)


def block_averaging(cutoff, block_size, data):
    data_cut = data[cutoff:]
    data_len = len(data_cut)
    number_of_blocks = int(data_len/block_size)
    average_array = np.zeros(number_of_blocks)
    for i in range(number_of_blocks):
        average_array[i] = (data_cut[i * block_size:(i + 1) * block_size]).mean()
    averages = average_array.mean()
    stdv = np.std(average_array, ddof=1) / np.sqrt(number_of_blocks)

    return number_of_blocks, averages, stdv


def initial_position(mass, w, hbar, N_particles, beads):
    mu, sigma = 0, math.sqrt(hbar / (mass * w))
    x = np.random.normal(mu, sigma, N_particles * beads)
    return x


def initial_velocity(mass, beta, N_particles, beads):
    mu, sigma = 0, math.sqrt(1 / (beta * mass))
    vx = np.random.normal(mu, sigma, N_particles * beads)
    return vx


def oscillator_potential(mass, w, x):
    potential = (0.5 * mass * w**2 * x**2).sum()
    return potential


def springs_potential(mass, wp, x):
    potential = 0
    for j in range(len(x)):
        if j == (len(x) - 1):
            potential += 0.5 * mass * wp ** 2 * (x[0] - x[j]) ** 2
        else:
            potential += 0.5 * mass * wp ** 2 * (x[j + 1] - x[j]) ** 2
    return potential


def oscillator_force(mass, w, x):
    f1 = - mass * w ** 2 * x
    return f1


def ring_springs_force(mass, wp, x):
    c2 = - mass * wp ** 2
    f2 = np.zeros(len(x))
    for j in range(len(x)):
        if j == (len(x) - 1):
            f2[j] = c2 * (2 * x[j] - x[0] - x[j - 1])
        elif j == 0:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[-1])
        else:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[j - 1])
    return f2


def kinetic_1d(mass, v):
    kinetic = (0.5 * mass * v ** 2).sum()
    return kinetic


def gamma(dt):  # This is the friction term
    g = 1 / (100 * dt)
    return g


def xsi(beads):  # This is the noise term
    z = np.random.normal(0, 1, beads)
    return z


def langevin(mass, beta, v, dt, beads, N_particles):
    kin1 = kinetic_1d(mass, v)  # Plus for Ethermo
    c1 = math.exp(-1 * gamma(dt) * dt / 2)
    c2 = math.sqrt((1 / beta / mass) * (1 - c1 ** 2))
    vel = c1 * v + c2 * xsi(beads * N_particles)
    kin2 = kinetic_1d(mass, vel)  # Minus for Ethermo
    return vel, kin1, kin2


# def kinetic_estimator_primitive(beads, beta, v_oo, v_o):
#     exp_v_oo = np.exp(- beta * v_oo)
#     exp_v_o = np.exp(- beta * v_o)
#     k_estimator = (beads / beta) - ((exp_v_oo * v_oo + exp_v_o * v_o) / (exp_v_oo + exp_v_o))
#     return k_estimator


def kinetic_estimator_virial2(beads, mass, w, x):
    k_estimator = - (1 / (2 * beads)) * (x * (- mass * w**2 * x)).sum()
    return k_estimator


def force_cumulative(mass, wp, beta, beads, N_particles, x_long):
    x = np.array(np.split(x_long, N_particles))
    # Potential
    v_o = springs_potential(mass, wp, x_long)  # One long ring
    v_oo = springs_potential(mass, wp, x[0]) + springs_potential(mass, wp, x[1])  # Two rings
    exp_v_o = np.exp(- beta * (v_o - v_oo))
    exp_v_oo = np.exp(- beta * (v_oo - v_oo))
    z_exp = exp_v_oo + exp_v_o
    # Forces
    f_oo_1 = (ring_springs_force(mass, wp, x[0])) * exp_v_oo
    f_oo_2 = (ring_springs_force(mass, wp, x[1])) * exp_v_oo
    f_o = (ring_springs_force(mass, wp, x_long)) * exp_v_o
    # Comulative Avg Force
    f1 = (f_o[:beads] + f_oo_1)
    f2 = (f_o[beads:] + f_oo_2)
    f = np.concatenate((f1, f2)) / z_exp
    return f, z_exp, v_o, v_oo


def langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, beads, N_particles):
    t = 0
    s = 0
    times = np.zeros(g_steps)
    steps = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    potential = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    e_change = np.zeros(g_steps)
    h_eff = np.zeros(g_steps)
    h_eff_change = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    pot_est = np.zeros(g_steps)
    kin_est = np.zeros(g_steps)
    # Histogram
    pos = np.zeros(g_steps)
    vel = np.zeros(g_steps)

    x = initial_position(mass, w, hbar, N_particles, beads)
    vx = initial_velocity(mass, beta, N_particles, beads)
    f_weighted, z_exp, v_o, v_oo = force_cumulative(mass, wp, beta, beads, N_particles, x)
    f_ext = oscillator_force(mass, w, x) / beads
    force = f_weighted + f_ext
    ethermo = 0.0
    for step in range(0, g_steps):
        times[step] = t
        steps[step] = int(s)
        kin[step] = kinetic_1d(mass, vx)
        potential[step] = v_oo - (1 / beta) * np.log((z_exp/2)) + oscillator_potential(mass, w, x) / beads
        e_tot[step] = kin[step] + potential[step]
        e_change[step] = abs((e_tot[step] - e_tot[0])) * 100 / e_tot[0]  # abs?
        h_eff[step] = e_tot[step] + ethermo
        h_eff_change[step] = abs((h_eff[step] - h_eff[1000]) * 100) / h_eff[1000]
        temp_exp[step] = 2.0 * kin[step] / (kboltz * beads * N_particles)  # Divided by beads from the Equipartition Function
        pot_est[step] = oscillator_potential(mass, w, x) / beads
        # kin_est[step] = kinetic_estimator_primitive(beads, beta, v_oo, v_o)
        kin_est[step] = kinetic_estimator_virial2(beads, mass, w, x)
        t += dt
        s += 1
        # Histogram
        # pos[step] = x[0]
        # vel[step] = vx[0]

        # Time propagation
        vx, kin1, kin2 = langevin(mass, beta, vx, dt, beads, N_particles)
        ethermo += kin1 - kin2
        vx = vx + 0.5 * dt * (force / mass)
        x = x + vx * dt
        f_weighted, z_exp, v_o, v_oo = force_cumulative(mass, wp, beta, beads, N_particles, x)
        f_ext = oscillator_force(mass, w, x) / beads
        force = f_weighted + f_ext
        vx = vx + 0.5 * dt * (force / mass)
        vx, kin1, kin2 = langevin(mass, beta, vx, dt, beads, N_particles)
        ethermo += kin1 - kin2
    return steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change


def force_cumulative_fermion(mass, wp, beta, beads, N_particles, x_long):
    x = np.array(np.split(x_long, N_particles))
    # Potential
    v_o = springs_potential(mass, wp, x_long)  # One long ring
    v_oo = springs_potential(mass, wp, x[0]) + springs_potential(mass, wp, x[1])  # Two rings
    exp_v_o = np.exp(- beta * (v_o - v_oo))
    exp_v_oo = np.exp(- beta * (v_oo - v_oo))
    z_exp = exp_v_oo + exp_v_o
    # Forces
    f_oo_1 = (ring_springs_force(mass, wp, x[0])) * exp_v_oo
    f_oo_2 = (ring_springs_force(mass, wp, x[1])) * exp_v_oo
    f_o = (ring_springs_force(mass, wp, x_long)) * exp_v_o
    # Comulative Avg Force
    f1 = (f_o[:beads] + f_oo_1)
    f2 = (f_o[beads:] + f_oo_2)
    f = np.concatenate((f1, f2)) / z_exp
    return f, z_exp, v_o, v_oo

def langevin_dynamics_fermion(g_steps, dt, mass, beta, hbar, kboltz, w, wp, beads, N_particles):
    t = 0
    s = 0
    times = np.zeros(g_steps)
    steps = np.zeros(g_steps)
    kin = np.zeros(g_steps)
    potential = np.zeros(g_steps)
    e_tot = np.zeros(g_steps)
    e_change = np.zeros(g_steps)
    h_eff = np.zeros(g_steps)
    h_eff_change = np.zeros(g_steps)
    temp_exp = np.zeros(g_steps)
    pot_est = np.zeros(g_steps)
    kin_est = np.zeros(g_steps)
    s_o_dinom = np.zeros(g_steps)
    s_o_num = np.zeros(g_steps)
    # Histogram
    pos = np.zeros(g_steps)
    vel = np.zeros(g_steps)

    x = initial_position(mass, w, hbar, N_particles, beads)
    vx = initial_velocity(mass, beta, N_particles, beads)
    f_weighted, z_exp, v_o, v_oo = force_cumulative(mass, wp, beta, beads, N_particles, x)
    f_ext = oscillator_force(mass, w, x) / beads
    force = f_weighted + f_ext
    ethermo = 0.0
    for step in range(0, g_steps):
        times[step] = t
        steps[step] = int(s)
        kin[step] = kinetic_1d(mass, vx)
        potential[step] = v_oo - (1 / beta) * np.log((z_exp/2)) + oscillator_potential(mass, w, x) / beads
        e_tot[step] = kin[step] + potential[step]
        e_change[step] = abs((e_tot[step] - e_tot[0])) * 100 / e_tot[0]  # abs?
        h_eff[step] = e_tot[step] + ethermo
        h_eff_change[step] = abs((h_eff[step] - h_eff[1000]) * 100) / h_eff[1000]
        temp_exp[step] = 2.0 * kin[step] / (kboltz * beads * N_particles)  # Divided by beads from the Equipartition Function
        pot_est[step] = oscillator_potential(mass, w, x) / beads
        # kin_est[step] = kinetic_estimator_primitive(beads, beta, v_oo, v_o)
        kin_est[step] = kinetic_estimator_virial2(beads, mass, w, x)
        s_o_dinom[step] = (np.exp(-beta * v_o) - np.exp(-beta * v_oo)) / (np.exp(-beta * v_o) + np.exp(-beta * v_oo))
        s_o_num[step] = (pot_est[step] * s_o_dinom[step])
        t += dt
        s += 1
        # Histogram
        # pos[step] = x[0]
        # vel[step] = vx[0]

        # Time propagation
        vx, kin1, kin2 = langevin(mass, beta, vx, dt, beads, N_particles)
        ethermo += kin1 - kin2
        vx = vx + 0.5 * dt * (force / mass)
        x = x + vx * dt
        f_weighted, z_exp, v_o, v_oo = force_cumulative(mass, wp, beta, beads, N_particles, x)
        f_ext = oscillator_force(mass, w, x) / beads
        force = f_weighted + f_ext
        vx = vx + 0.5 * dt * (force / mass)
        vx, kin1, kin2 = langevin(mass, beta, vx, dt, beads, N_particles)
        ethermo += kin1 - kin2
    return steps, times, pos, vel, kin, potential, e_tot, e_change, \
           temp_exp, pot_est, kin_est, h_eff_change, s_o_dinom, s_o_num




