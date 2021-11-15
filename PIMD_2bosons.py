from MD_QHO_Functions_2bosons import *



# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, g, si, beads, N_particles)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=1000, data=temp_exp)
# number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=1000, data=pot_est)
# number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=1000, data=kin_est)
#
# # number_of_blocks, avg_etot, stdv_etot = block_averaging(cutoff, block_size=1000, data=e_tot)
# # number_of_blocks1, avg_potential_est, stdv_potential = block_averaging(cutoff, block_size=1000, data=potential)
# # number_of_blocks2, avg_kin_est, stdv_kin = block_averaging(cutoff, block_size=1000, data=kin)
#
# print("Mean Kinetic Estimator:", avg_kin_est, "+-", stdv_kin)
# print("Mean Potential Estimator:", avg_potential_est, "+-", stdv_potential)
# print("Mean Temperature:", avg_temp_exp, "+-", stdv_temp)
# print("Mean Total Energy Estimator :", avg_potential_est + avg_kin_est, "+-", np.sqrt(stdv_potential**2 + stdv_kin**2))
#
# np.savez("2B_QHO_beta10.1", avg_potential_est=avg_potential_est,
#          avg_kin_est=avg_kin_est, stdv_kin=stdv_kin, stdv_potential=stdv_potential)


beads_array = np.array([3, 8, 16, 24, 28, 32])

e_tot_array = np.zeros(len(beads_array))
e_tot_stdv_array = np.zeros(len(beads_array))

for i, b in enumerate(beads_array):
    wp = math.sqrt(b) / (beta * hbar)
    start = time.time()
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
        langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, g, si, b, N_particles)
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

print("e array:", e_tot_array)

np.savez("2B_QHO_C_beta10_vs_beads_20K",  e_tot_array=e_tot_array, e_tot_stdv_array=e_tot_stdv_array, beads_array=beads_array)

pass

figtotenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, e_tot_array, '.', label="Mean Total Energy", color="blue")
plt.errorbar(beads_array, e_tot_array, yerr=e_tot_stdv_array, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
plt.show()

numb = int(g_steps*0.3)*0
# figenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times[numb:], e_tot[numb:], '.', label="Total Energy vs Time", color="black")
# plt.plot(times[numb:], kin[numb:], '.', label="Kinetic vs Time", color="red")
# plt.plot(times[numb:], potential[numb:], '.', label="Potential vs Time", color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()
#
# figper = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times[numb:], e_change[numb:], '.', label="Energy % Change", color="red")
# plt.xlabel("time")
# plt.ylabel("(E(t) - E0)*100 / E0  [%]")
# # plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()

figtemp = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times[numb:], temp_exp[numb:], '.', label="Temperature", color="red")
plt.xlabel("time")
plt.ylabel("Temperature")
plt.legend()
plt.show()

figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times[numb:], h_eff_change[numb:], '.', label="Effective Hamiltonian", color="red")
plt.xlabel("time")
plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
# plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()

figpotest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(steps[numb:], pot_est[numb:], '.', label="P_Estimator vs step", color="blue")
plt.xlabel("Steps")
plt.ylabel("Potential Estimator")
plt.legend()
plt.show()
#
figkinest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(steps[numb:], kin_est[numb:], '.', label="K_Estimator vs step", color="blue")
plt.xlabel("Steps")
plt.ylabel("Kinetic Estimator")
plt.legend()
plt.show()



 # Results


# betas = np.array([1,  2,  3,  5, 10])
# e_tot_b = np.array([1.98825745, 1.23819174, 1.08257483, 1.01779556, 0.99359082])
# stdv_b = np.array([0.11363226, 0.05506471, 0.03580274, 0.02237071, 0.01475264])

# beta= 2 60K itr 1.190640207206084 +- 0.03331297434928319
# beta = 2 20K itr 1.23819174 +- 0.05506471
# beta = 0.5   20K itr  3.6827085325 +- 0.22941947274798993
# beta = 0.5   60K itr 3.71512931 +- 0.12919464121210192
betas = np.array([0.5, 1,  2,  3,  5, 10])
e_tot_b = np.array([3.71512931, 1.98825745, 1.190640207206084, 1.08257483, 1.01779556, 0.99359082])
stdv_b = np.array([0.12919464121210192, 0.11363226, 0.03331297434928319, 0.03580274, 0.02237071, 0.01475264])

# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# q = np.linspace(0.3, 11, 1000)
# p = (hbar * w * (np.exp(q * hbar * w) + np.exp(2*q * hbar * w) + 2))/(np.exp(2 * q * hbar * w)-1)
# plt.plot(q, p, 'g',  label="Analytical Result")
# plt.plot(betas, e_tot_b, '.', label="Mean Total Energy", color="blue", linestyle='None')
# plt.errorbar(betas, e_tot_b, yerr=stdv_b, ecolor="black", linestyle='None')
# plt.xlabel("Beta")
# plt.ylabel("Mean Total Energy")
# plt.legend(loc='upper right')
# plt.show()


# e05 =[3.51364408, 3.53553696, 3.74749349, 3.78348673, 3.42482214, 3.47558984, 3.70211219, 3.67790933, 3.65613724, 3.69467537]
# e1 =[1.77756524, 1.88358542, 1.97333656, 2.02469105, 1.81690049, 2.05772718, 2.00220635, 1.96750464, 1.99506137]
# e2 =[1.0917804,  1.17912541, 1.23402858, 1.251729, 1.14492413, 1.24778852, 1.24208449, 1.23004114, 1.2424496]
# e3 = [0.91550385, 1.03109412, 1.07090572, 1.09091144, 1.01400037, 1.08825822, 1.07807176, 1.08580715, 1.08384558]
# e5 = [0.74800101, 0.93259973, 0.99202508, 1.01461955, 0.96732838, 1.01736942, 1.0201008, 1.01360915, 1.01967671]
# e10 = [0.49658288, 0.77110258, 0.90809851, 0.96781267, 0.94231154, 0.97672965, 1.00596272, 0.98165858, 0.99315118]
# s05 = [0.04742002, 0.05975375, 0.07407656, 0.0766224, 0.08601928, 0.10037365, 0.1072404, 0.11472529, 0.1143131, 0.12207954]
# s1 =[0.02245057, 0.03421171, 0.03982982, 0.04471548, 0.04997226, 0.0642453, 0.0592174, 0.06288521, 0.07383116]
# s2 = [0.01046913, 0.01581804, 0.01978168, 0.02154474, 0.02222855, 0.03020612, 0.02983664, 0.02898888, 0.03607689]
# s3 = [0.00710411, 0.01090228, 0.01325285, 0.01464614, 0.01512544, 0.02061812, 0.01927779, 0.01948213, 0.02303583]
# s5 = [0.00500744, 0.00720508, 0.00929217, 0.0102218, 0.01001999, 0.01265714, 0.01269869, 0.0120366, 0.01393959]
# s10 = [0.00331076, 0.00521686, 0.00655702, 0.00705895, 0.0067865, 0.00805114, 0.00880198, 0.00765436, 0.00903196]

# beads_array = np.array([3, 5, 8, 10, 16, 24, 28, 32, 36])
# e1 =[1.8833593120295147, 1.9245363730641891, 1.8838085992454328, 1.9019961437915973, 1.959432763979145, 1.9095407697127011, **2.00220635, **1.96750464, **1.99506137]
# s1 =[0.011009353908803352, 0.013850426694430757, 0.015995631073664696, 0.017808787453056422, 0.022128074438208238,  0.024033140264882163, **0.0592174, **0.06288521, **0.07383116]