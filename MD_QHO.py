import matplotlib.pyplot as plt
import numpy as np

from MD_QHO_Functions import *

# CONSTANTS
mass = 1
w = 1
hbar = 1
kboltz = 1
beta = 10
beta_array = np.array([1, 2, 3, 6, 8, 10])

beads = 32
beads_array = np.array([2, 8, 16, 32])

# Time
T_max = 2000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)

# ErrorBars
cutoff = int(g_steps * 0.3)
num_blocks = 3

seed = 347885
np.random.seed(seed)

# Main

# Calculated Mean Estimators (kinetic or potential) for a given number of beads

# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)
# # avg, stdv = block_averaging(g_steps, num_of_blocks=3, data=temp_exp)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
# print("Mean Temperature:", np.mean(temp_exp[cutoff:]))
# print("Mean Potential_estimator:", np.mean(pot_est[cutoff:]))
# print("Mean Kinetic_estimator:", np.mean(kin_est[cutoff:]))


# See the Fluctuations in the Potential Estimator and Kinetic Estimator (Primitive) vs time
# Different number of beads

# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est0, kin_est0, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array[0])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est1, kin_est1, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array[1])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est2, kin_est2, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array[2])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est3, kin_est3, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array[3])
# stop = time.time()
# duration = stop - start


# Different Temperatures - <Etot_est> vs time

# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est0, kin_est0, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[0], hbar, kboltz, w, beads)
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est1, kin_est1, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[1], hbar, kboltz, w, beads)
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est2, kin_est2, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[2], hbar, kboltz, w,beads)
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est3, kin_est3, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[3], hbar, kboltz, w,beads)
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est4, kin_est4, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[4], hbar, kboltz, w, beads)
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est5, kin_est5, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta_array[5], hbar, kboltz, w, beads)
#
# stop = time.time()
# e_tot_est0 = np.mean(kin_est0[cutoff:]) + np.mean(pot_est0[cutoff:])
# e_tot_est1 = np.mean(kin_est1[cutoff:]) + np.mean(pot_est1[cutoff:])
# e_tot_est2 = np.mean(kin_est2[cutoff:]) + np.mean(pot_est2[cutoff:])
# e_tot_est3 = np.mean(kin_est3[cutoff:]) + np.mean(pot_est3[cutoff:])
# e_tot_est4 = np.mean(kin_est4[cutoff:]) + np.mean(pot_est4[cutoff:])
# e_tot_est5 = np.mean(kin_est5[cutoff:]) + np.mean(pot_est5[cutoff:])
#
# e_tot_est_array = np.array([e_tot_est0, e_tot_est1, e_tot_est2, e_tot_est3, e_tot_est4, e_tot_est5])
# duration = stop - start
# print("Duration time:", duration)


# <Etot_est> vs Beads>

# start = time.time()
# mean_e_tot_est, stdv = \
#     langevin_dynamics_beads(g_steps, cutoff, num_blocks, dt, mass, beta, hbar, kboltz, w, beads_array)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
# print("Mean ", mean_e_tot_est)
# print("STDV", stdv)




# Plots

# Histograms
# plt.figure(1)
# n, bins, patches = plt.hist(pos, 50, density=True, facecolor='b', alpha=0.75)
# plt.xlabel('Position')
# plt.ylabel('Counts')
# plt.title('Position of Bead')
# plt.grid(True)
# plt.show()

# plt.figure(2)
# n, bins, patches = plt.hist(vel, 50, density=True, facecolor='b', alpha=0.75)
# plt.xlabel('Velocity')
# plt.ylabel('Counts')
# plt.title('Velocity of Bead')
# plt.grid(True)
# plt.show()


# Conservation of Energy and Temperature
# print("Mean Total Energy:", np.mean(e_tot[cutoff:]))

# figtemp = plt.figure()
# plt.plot(times, temp_exp, '.', label="Temperature" , color="red")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()

# figper = plt.figure()
# plt.plot(times, e_change, '.', label="Energy % Change" , color="red")
# plt.xlabel("time")
# plt.ylabel("(E(t) - E0)*100 / E0  [%]")
# plt.ylim([-0.1, 0. 1])
# plt.legend()
# plt.show()

# figenergy = plt.figure()
# plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
# plt.plot(times, potential, '.', label="Potential vs Time" , color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()

# Estimators

# figpotest = plt.figure()
# plt.plot(steps, pot_est, '.', label="P_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()

# figenest = plt.figure()
# plt.plot(beads_array, mean_e_tot_est, '.', label="Mean Total Energy Estimator", color="blue")
# plt.errorbar(beads_array, mean_e_tot_est, yerr=stdv, ecolor="black")
# plt.xlabel("Beads")
# plt.ylabel("Mean Tot_Energy Estimator")
# plt.legend()
# plt.show()

# figkinest = plt.figure()
# plt.plot(times, kin_est0, '.', label="Kin_Est 0", color="black")
# plt.plot(times, kin_est1, '.', label="Kin_Est 1", color="red")
# plt.plot(times, kin_est2, '.', label="Kin_Est 2", color="blue")
# plt.plot(times, kin_est3, '.', label="Kin_Est 3", color="green")
# plt.xlabel("time")
# plt.ylabel("Kinetic Energy Estimator")
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.plot(times, pot_est0, '.', label="Pot_Est 0", color="black")
# plt.plot(times, pot_est1, '.', label="Pot_Est 1", color="red")
# plt.plot(times, pot_est2, '.', label="Pot_Est 2", color="blue")
# plt.plot(times, pot_est3, '.', label="Pot_Est 3", color="green")
# plt.xlabel("time")
# plt.ylabel("Potential Energy Estimator")
# plt.legend()
# plt.show()

figetotest = plt.figure()
plt.plot(beta_array, e_tot_est_array, label="Mean Total Energy Estimator", color="black")
plt.xlabel("beta")
plt.ylabel("Mean Total Energy Estimator")
plt.legend()
plt.show()


# Effective Hemiltonian for Energy Conservation

# figper = plt.figure()
# plt.plot(times, h_eff_change, '.', label="Effective Hemiltonian", color="red")
# plt.xlabel("time")
# plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
# plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()


