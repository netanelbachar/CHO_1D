from MD_QHO_Functions_2bosons import *


# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, beads, N_particles)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=1000, data=temp_exp)
# number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=3000, data=pot_est)
# number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=3000, data=kin_est)
#
# print("Mean Temperature:", avg_temp_exp, "+-", stdv_temp)
# print("Mean Potential_estimator:", avg_potential_est, "+-", stdv_potential)
# print("Mean Kinetic_estimator:", avg_kin_est, "+-", stdv_kin)
# print("Mean Total Energy Estimator:", avg_potential_est + avg_kin_est, "+-", np.sqrt(stdv_potential**2 + stdv_kin**2))
#
# np.savez("2B_QHO_beta10.1", avg_potential_est=avg_potential_est,
#          avg_kin_est=avg_kin_est, stdv_kin=stdv_kin, stdv_potential=stdv_potential)


beads_array = np.array([3, 5, 8, 10, 16, 20, 24, 28, 32, 36])
# beads_array = np.array([3, 6, 10, 16, 20, 24, 28, 32, 36])
e_tot_array = np.zeros(len(beads_array))
e_tot_stdv_array = np.zeros(len(beads_array))
for i, b in enumerate(beads_array):
    wp = math.sqrt(b) / (beta * hbar)

    start = time.time()
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
        langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, b, N_particles)
    stop = time.time()
    duration = stop - start
    print("Time of execution:", duration)
    print("Iteration:", i)
    number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=2000, data=temp_exp)
    number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=6000, data=pot_est)
    number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=6000, data=kin_est)

    print("Mean classical Temperature:", avg_temp_exp, "+-", stdv_temp)
    print("Mean Potential_estimator:", avg_potential_est, "+-", stdv_potential)
    print("Mean Kinetic_estimator:", avg_kin_est, "+-", stdv_kin)
    print("Mean Total Energy Estimator:", avg_potential_est + avg_kin_est, "+-",
          np.sqrt(stdv_potential ** 2 + stdv_kin ** 2))

    e_tot_array[i] = avg_potential_est + avg_kin_est
    e_tot_stdv_array[i] = np.sqrt(stdv_potential ** 2 + stdv_kin ** 2)

print ("e array:", e_tot_array)

np.savez("2B_QHO_beta0.5_vs_beads_20K",  e_tot_array=e_tot_array, e_tot_stdv_array=e_tot_stdv_array, beads_array=beads_array)

pass

figtotenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, e_tot_array, '.', label="Mean Total Energy", color="blue")
plt.errorbar(beads_array, e_tot_array, yerr=e_tot_stdv_array, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
plt.show()



figenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, e_tot, '.', label="Total Energy vs Time", color="black")
plt.plot(times, kin, '.', label="Kinetic vs Time", color="red")
plt.plot(times, potential, '.', label="Potential vs Time", color="blue")
plt.xlabel("time")
plt.ylabel("Energy")
plt.legend()
plt.show()

figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, e_change, '.', label="Energy % Change", color="red")
plt.xlabel("time")
plt.ylabel("(E(t) - E0)*100 / E0  [%]")
plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()



figtemp = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, temp_exp, '.', label="Temperature", color="red")
plt.xlabel("time")
plt.ylabel("Temperature")
plt.legend()
plt.show()
#
figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, h_eff_change, '.', label="Effective Hamiltonian", color="red")
plt.xlabel("time")
plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()

figpotest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(steps, pot_est, '.', label="P_Estimator vs step", color="blue")
plt.xlabel("Steps")
plt.ylabel("Potential Estimator")
plt.legend()
plt.show()

figkinest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
plt.xlabel("Steps")
plt.ylabel("Kinetic Estimator")
plt.legend()
plt.show()


# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# q = np.linspace(0.8, 11, 1000)
# p = ((2*hbar * w * np.exp(q * hbar * w))*(np.exp(q * hbar * w)+np.exp(2* q * hbar * w)+1))/((np.exp(q * hbar * w)-1)*(np.exp(q * hbar * w)+1)*(2*np.exp(q * hbar * w)+1))
# plt.plot(q, p, 'g')
# plt.plot(beta, e_tot_b, '.', label="Mean Total Energy", color="blue", linestyle='None')
# plt.errorbar(beta, e_tot_b, yerr=stdv_b, ecolor="black", linestyle='None')
# plt.xlabel("Beta")
# plt.ylabel("Mean Total Energy")
# plt.legend(loc='upper right')
# plt.show()

# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# q = np.linspace(0.8, 11, 1000)
# p = hbar * w * (1 + 2 / (np.exp(q * hbar * w) - 1))
# plt.plot(q, p, 'g')
# plt.plot(beta, e_tot_b, '.', label="Mean Total Energy", color="blue", linestyle='None')
# plt.errorbar(beta, e_tot_b, yerr=stdv_b, ecolor="black", linestyle='None')
# plt.xlabel("Beta")
# plt.ylabel("Mean Total Energy")
# plt.legend(loc='upper right')
# plt.show()
