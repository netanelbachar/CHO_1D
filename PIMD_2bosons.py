from MD_QHO_Functions_2bosons import *


# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, wp_long, beads, N_particles)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
# number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=1000, data=temp_exp)
# number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=4000, data=pot_est)
# number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=4000, data=kin_est)

beads_array = np.array([55])
e_tot_array = np.zeros(len(beads_array))
e_tot_stdv_array = np.zeros(len(beads_array))
for i, b in enumerate(beads_array):
    wp = math.sqrt(b) / (beta * hbar)
    wp_long = math.sqrt(b * N_particles) / (beta * hbar)

    start = time.time()
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
        langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, wp_long, b, N_particles)
    stop = time.time()
    duration = stop - start
    print("Time of execution:", duration)

    number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=1000, data=temp_exp)
    number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=4000, data=pot_est)
    number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=4000, data=kin_est)

    print("Mean classical Temperature:", avg_temp_exp, "+-", stdv_temp)
    print("Mean Potential_estimator:", avg_potential_est, "+-", stdv_potential)
    print("Mean Kinetic_estimator:", avg_kin_est, "+-", stdv_kin)
    print("Mean Total Energy Estimator:", avg_potential_est + avg_kin_est)

    e_tot_array[i] = avg_kin_est + avg_potential_est
    e_tot_stdv_array[i] = np.sqrt(stdv_potential**2 + stdv_kin**2)

print ("e array:", e_tot_array)

figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, e_tot_array, '.', label="total energy", color="red")
plt.xlabel("beads")
plt.ylabel("Total Energy")
plt.legend()
plt.show()

# figenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, e_tot, '.', label="Total Energy vs Time", color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time", color="red")
# plt.plot(times, potential, '.', label="Potential vs Time", color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()
#
# figper = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, e_change, '.', label="Energy % Change", color="red")
# plt.xlabel("time")
# plt.ylabel("(E(t) - E0)*100 / E0  [%]")
# plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()


# figtemp = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, temp_exp, '.', label="Temperature", color="red")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()

# figper = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, h_eff_change, '.', label="Effective Hamiltonian", color="red")
# plt.xlabel("time")
# plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
# plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(steps, pot_est, '.', label="P_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()



