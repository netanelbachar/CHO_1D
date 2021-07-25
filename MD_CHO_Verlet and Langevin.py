# Classical Harmonic Oscilator in One Dimension
from functions_QHO import *

# Constants
mass = 1
w = 10
kboltz = 1
hbar = 1
temperature = 1
beta = 1 / (temperature * kboltz)
beads = 1  # Only 1 Particle
seed = 347885
np.random.seed(seed)

# Time
T_max = 100
dt = 1 * 10 ** (-4)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)


# Main Code
algorithm = input("v for verlet or l from langevin:")
if algorithm == "v":

    times, pos, vx_n, kin, pot, e_tot, temp_exp = verlet_classic(g_steps, mass, dt, beta, w, hbar, kboltz, beads)

    plt.figure(figsize=[10, 8])
    n, bins, patches = plt.hist(x=pos, bins=30, density=True, alpha=0.7,rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Position', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.title('MD_CHO NVE : Position Histogram', fontsize=15)
    plt.show()

    fig3 = plt.figure()
    plt.plot(times, e_tot, label="Total Energy", color="red")
    plt.plot(times, pot, label="Pot", color="blue")
    plt.plot(times, kin, label="kin", color="black")
    plt.ylim(0, 100)
    plt.xlabel("t")
    plt.ylabel("Energy")
    plt.legend()
    plt.show()

else:
    times, pos, vx_n, kin, pot, e_tot, temp_exp = langevin_classic(g_steps, mass, dt, beta, w, hbar, kboltz, beads)
    print("Mean Temperature:", sum(temp_exp) / len(temp_exp))
# PLOT NVT
    plt.figure(figsize=[10, 8])
    n, bins, patches = plt.hist(x=pos, bins=50, density=True, alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Position',fontsize=15)
    plt.ylabel('Count',fontsize=15)
    plt.title('MD_CHO NVT : Position Histogram',fontsize=15)
    plt.show()

    plt.figure(figsize=[10, 8])
    n, bins, patches = plt.hist(x=vx_n, bins=50, density=True, alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Velocity', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.title('MD_CHO NVT : Velocity Histogram', fontsize=15)
    plt.show()

#     fig2 = plt.figure()
#     plt.plot(times, temp_exp, label="Temperature", color="red")
#     plt.xlabel("Time")
#     plt.ylabel("Temperature")
#     plt.legend()
#     plt.show()
#
#
# fig2 = plt.figure()
# plt.plot(times, pos, label="Position", color="red")
# plt.xlabel("Time")
# plt.ylabel("Position")
# plt.legend()
# plt.show()











# ###############NOT TOUCH
#
# # Quantum Harmonic Oscillator in One Dimension
#
# # IMPORTS
# import numpy as np
# import scipy as sp
# import matplotlib.pyplot as plt
# import math
# # from functions_QHO import *
#
# # CONSTANTS
# N_particles = 1
# mass = 1
# w = 10
# kboltz = 1  # 1.38064852 * 10 ** (-23) # [J/K]
# hbar = 1  # 6.62607004 * 10 ** (-34)  # [J s] Planck Constant
# temperature = 10  # 298.15  # [K]
# beta = 1 / (temperature * kboltz)
# beads = 8
# wp = math.sqrt(beads) / (beta * hbar)
# # Time
# T_max = 1
# dt = 1 * 10 ** (-4)
# g_steps = int(T_max / dt)
# print("Number of Steps:", g_steps)
# seed = 3477
# np.random.seed(seed)
#
# def initial_position():
#     """
#
#     :return: Position Array of zeros with the first bead at different position
#     length as the number of beads
#     """
#
#     x = np.zeros(beads)
#     x[0] = 1
#     return x
#
#
# def velocity():
#     """
#
#     :return: Velocity Array of zeros.
#     length as the number of beads
#     """
#     vx = np.zeros(beads)
#     return vx
#
#
# def force(x):
#     """
#
#     :param x: array of bead positions in 1D
#     :return: array of forces acting on each bead
#     """
#     force_beads = np.zeros(len(x))
#     for j in range(len(x)):
#         if j == (len(x) - 1):
#             force_beads[j] = - mass * wp**2 * (2 * x[j] - x[0] - x[j-1]) - (mass * w**2 / beads) * x[j]
#         elif j == 0:
#             force_beads[j] = - mass * wp**2 * (2 * x[j] - x[j+1] - x[-1]) - (mass * w**2 / beads) * x[j]
#         else:
#             force_beads[j] = - mass * wp**2 * (2 * x[j] - x[j + 1] - x[j-1]) - (mass * w**2 / beads) * x[j]
#     return force_beads
#
#
# def potential_1d(x):
#     """
#
#     :param x: Array of bead positions
#     :return: Sum of all potential energy of the beads
#     """
#     potential = 0
#     for j in range(len(x)):
#         if j == (len(x) - 1):
#             potential += 0.5 * mass * wp**2 * (x[0] - x[j])**2 + ((mass * w**2) / (2 * beads)) * x[j]**2
#         else:
#             potential += 0.5 * mass * wp**2 * (x[j+1] - x[j])**2 + ((mass * w**2) / (2 * beads)) * x[j]**2
#     return potential
#
#
# def kinetic_1d(v):
#     """
#
#     :param v: Velocity array of the beads
#     :return: Total Kinetic energy of the beads
#     """
#     kinetic = 0
#     for j in range(len(v)):
#         kinetic += 0.5 * mass * v[j]**2
#     return kinetic
#
#
# # Langevin Dynamics
# def l_gamma(dt):  # This is the friction term
#     gamma = 0.01 / dt
#     return gamma
#
#
# def l_xsi():  # This is the noise term
#     xsi = np.random.normal(0, 1, 1)
#     return xsi
#
#
# def langevin(v, time_step):
#     vel = math.exp(-1 * l_gamma(time_step) * time_step / 2) * v + \
#           (1 / (math.sqrt(mass * beta))) * math.sqrt(1-math.exp(-1 * l_gamma(time_step) * time_step)) * l_xsi()
#     return vel
#
#
# # Main Code
# t = 0
# x = initial_position()
# vx = velocity()
# force_x = force(x)
#
# times = np.zeros(g_steps)
# pos = np.zeros(g_steps)
# pos1 = np.zeros(g_steps)
# pos2 = np.zeros(g_steps)
# pos3 = np.zeros(g_steps)
# pos4 = np.zeros(g_steps)
# pos5 = np.zeros(g_steps)
# vx_n = np.zeros(g_steps)
# kin = np.zeros(g_steps)
# pot = np.zeros(g_steps)
# e_tot = np.zeros(g_steps)
# percent_change = np.zeros(g_steps)
# temp_exp = np.zeros(g_steps)
#
#
# algorithm = input("v for verlet or l from langevin:")
# if algorithm == "v":
#
#     for step in range(0, g_steps):
#         times[step] = t
#         t += dt
#         kin[step] = kinetic_1d(vx)
#         pot[step] = potential_1d(x)
#         e_tot[step] = kin[step] + pot[step]
#         pos[step] = x[0]
#         percent_change[step] = abs(e_tot[step] - e_tot[0]) * 100 / e_tot[0]
#         # vx_n[step] = vx
#         temp_exp[step] = kboltz * e_tot[step] / beads
#         for i in range(len(x)):
#         # Verlet
#             vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
#             x[i] = x[i] + dt * vx[i]
#             force_x = force(x)
#             vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
# # HIST NVE
# #     plt.figure(figsize=[10, 8])
# #     n, bins, patches = plt.hist(x=pos, bins=30, density=True, alpha=0.7,rwidth=0.85)
# #     plt.grid(axis='y', alpha=0.75)
# #     plt.xlabel('Position', fontsize=15)
# #     plt.ylabel('Count', fontsize=15)
# #     plt.title('MD_CHO NVE : Position Histogram', fontsize=15)
# #     plt.show()
#
# else:
#     for step in range(0, g_steps):
#         times[step] = t
#         t += dt
#         kin[step] = kinetic_1d(vx)
#         pot[step] = potential_1d(x)
#         e_tot[step] = kin[step] + pot[step]
#         # pos1[step] = x[0]
#         # pos2[step] = x[1]
#         # pos3[step] = x[2]
#         # pos4[step] = x[3]
#         # pos5[step] = x[4]
#         vx_n[step] = vx[0]
#         percent_change[step] = abs(e_tot[step] - e_tot[0]) * 100 / e_tot[0]
#         temp_exp[step] = kboltz * e_tot[step] / beads  # Divided by beads from the Equipartition Function
#
#         # Langevin
#         for i in range(len(x)):
#             vx[i] = langevin(vx[i], dt)
#             vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
#             x[i] = x[i] + dt * vx[i]
#             force_x = force(x)
#             vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
#             vx[i] = langevin(vx[i], dt)
#
# # PLOT NVT
#     plt.figure(figsize=[10, 8])
#     n, bins, patches = plt.hist(x=pos, bins=30, density=True, alpha=0.7, rwidth=0.85)
#     plt.grid(axis='y', alpha=0.75)
#     plt.xlabel('Position',fontsize=15)
#     plt.ylabel('Count',fontsize=15)
#     plt.title('MD_CHO NVT : Position Histogram',fontsize=15)
#     plt.show()
#
#
#     plt.figure(figsize=[10, 8])
#     n, bins, patches = plt.hist(x=vx_n, bins=30, density=True, alpha=0.7, rwidth=0.85)
#     plt.grid(axis='y', alpha=0.75)
#     plt.xlabel('Velocity', fontsize=15)
#     plt.ylabel('Count', fontsize=15)
#     plt.title('MD_CHO NVT : Velocity Histogram', fontsize=15)
#     plt.show()
#
# num = int(g_steps * 0.1)
# print("Mean Temperature:", sum(temp_exp[num:])/len(temp_exp[num:]))
#
#
#
# # PLOTS
#
#
# figper = plt.figure()
# plt.plot(times[:], percent_change[:], label="Percent Change", color="red")
# plt.xlabel("dt")
# plt.ylabel("% Change ((Ei-E0)*100/E0) ")
# # plt.ylim([-0.5, 0.5])
# plt.legend()
# plt.show()
#
#
# figtemp = plt.figure()
# plt.plot(times[:], temp_exp[:], '.', label="Temperature vs time" , color="black")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()
#
# # Position of particle in the ring
# # fig1 = plt.figure()
# # plt.plot(times, pos, '.', label="Position of Particle" , color="black")
# # plt.xlabel("time")
# # plt.legend()
# # plt.show()
#
# fig2 = plt.figure()
# plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
# plt.plot(times, pot, '.', label="Potential vs Time" , color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()
#
# # fig2 = plt.figure()
# # plt.plot(times, pos1, '.', label="Position 1" , color="black")
# # plt.plot(times, pos2, '.', label="Position 2" , color="red")
# # plt.plot(times, pos3, '.', label="Position 3", color="blue")
# # plt.plot(times, pos4, '.', label="Position 4", color="yellow")
# # plt.plot(times, pos5, '.', label="Position 5" , color="green")
# # plt.xlabel("time")
# # plt.ylabel("Position")
# # plt.legend()
# # plt.show()



##

#     for i in range(len(x)):
#         vx[i] = langevin(mass, beta, vx[i], dt)
#         vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
#         x[i] = x[i] + dt * vx[i]
#         force_x = force(mass, w, wp, beads, x)
#         vx[i] = vx[i] + 0.5 * dt * force_x[i] / mass
#         vx[i] = langevin(mass, beta, vx[i], dt)
# return times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp, pot_est

