import numpy as np

from MD_QHO_Functions import *

# Different Temperatures - <Etot_est> vs time

# What is left: Add STDV to the Graph

# bead_array_n = np.array([6, 8, 10, 14, 14, 20])
bead_array_n = np.array([5, 10, 10, 16, 20, 25])

start = time.time()

e_kin_est_array = np.zeros(len(bead_array_n))
e_pot_est_array = np.zeros(len(bead_array_n))
e_kin_est_stdv = np.zeros(len(bead_array_n))
e_pot_est_stdv = np.zeros(len(bead_array_n))

e_tot_est_array = np.zeros(len(bead_array_n))
e_tot_est_stdv = np.zeros(len(bead_array_n))

for i in range(0, len(bead_array_n)):
    print(beta_array[i])
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
        langevin_dynamics(g_steps, dt, mass, beta_array[i], hbar, kboltz, w, bead_array_n[i])
    e_kin_est_array[i], e_kin_est_stdv[i] = block_averaging(cutoff, num_blocks, data=kin_est)
    e_pot_est_array[i], e_pot_est_stdv[i] = block_averaging(cutoff, num_blocks, data=pot_est)
    e_tot_est_array[i], e_tot_est_stdv[i] = e_kin_est_array[i] + e_pot_est_array[i], e_kin_est_stdv[i] + e_pot_est_stdv[i]


stop = time.time()
print("Mean Tot Est:", e_tot_est_array, "+-", e_tot_est_stdv)
duration = stop - start
print("Duration time:", duration)

figetotest = plt.figure()
plt.rcParams.update({'font.size': 13})
q = np.linspace(0.8, 11, 1000)
p = hbar * w * (0.5 + 1 / (np.exp(q * hbar * w) - 1))
plt.rcParams.update({'font.size': 13})
plt.plot(beta_array, e_tot_est_array, label="Mean Total Energy Estimator", color="black")
plt.errorbar(beta_array, e_tot_est_array, yerr=e_tot_est_stdv, ecolor="black")
plt.plot(q, p, 'g')
plt.xlabel("beta")
plt.ylabel("Mean Total Energy Estimator")
plt.legend()
plt.show()


# bead_array_n = np.array([6, 8, 10, 14, 14, 20])
# Mean Tot Est: [0.97225615 0.59758398 0.63975371 0.50127964 0.50012068 0.48594599] +- [0.0842021  0.06610652 0.08265349 0.03628287 0.02547104 0.02999633]
# Duration time: 4021.7637960910797
#bead_array_n = np.array([5, 10, 10, 16, 20, 25])
# Mean Tot Est: [1.04053167 0.67612901 0.71810781 0.49913517 0.44919944 0.45556148] +- [0.12258558 0.07428768 0.08708906 0.05468467 0.03969176 0.03532206]
# Duration time: 4362.5687391757965


