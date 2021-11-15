from MD_QHO_Functions_2bosons import *

# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change, s_o_dinom, s_o_num, hminush, kin_est1 = \
#     langevin_dynamics_fermion_gaus(g_steps, dt, mass, beta, hbar, kboltz, w, wp, g, si, beads, N_particles)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=2000, data=temp_exp)
# number_of_blocks, kin_avg, stdv_kin = block_averaging(cutoff, block_size=4000, data=kin_est)
# number_of_blocks, potential_avg, stdv_potential = block_averaging(cutoff, block_size=4000, data=pot_est)
# number_of_blocks, avg_s_o_dinom, stdv_s_o_dinom = block_averaging(cutoff, block_size=4000, data=s_o_dinom)
# number_of_blocks, avg_s_o_num, stdv_s_o_num = block_averaging(cutoff, block_size=4000, data=s_o_num)
# number_of_blocks, avg_hminush, stdv_hminush = block_averaging(cutoff, block_size=4000, data=hminush)
#
# observable = (avg_s_o_num / avg_s_o_dinom)
#
# d_observable = (avg_s_o_num / avg_s_o_dinom) * \
#                np.sqrt((stdv_s_o_num / avg_s_o_num) ** 2 + (stdv_s_o_dinom / avg_s_o_dinom) ** 2)
#
# # E_H = E_H' + F-F'
# ferm_energy = observable + avg_hminush
#
# d_ferm_energy = math.sqrt(d_observable ** 2 + stdv_hminush ** 2)
#
# print("Mean Temperature:", avg_temp_exp, "+-", stdv_temp)
# print("Fermion: Total Energy Estimator:", observable, "+-", d_observable)
# print("<H - H'>_H' :", avg_hminush, "+-", stdv_hminush)
# print("Fermionic Energy :", ferm_energy, "+-", d_ferm_energy)



#
# # np.savez("2F_QHO_ensemble_beta05", observable=2*observable, d_observable=2*d_observable)
#
# plt.figure(1)
# plt.rcParams.update({'font.size': 13})
# n, bins, patches = plt.hist(pos, 100, density=True, facecolor='b', alpha=0.75)
# plt.xlabel('Position')
# plt.ylabel('Counts')
# plt.title('Position of Bead')
# plt.xlim([-10, 10])
# plt.grid(True)
# plt.show()

# figtemp = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, temp_exp, '.', label="Temperature", color="red")
# plt.xlabel("time")
# plt.ylabel("Temperature")
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
# figkinest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()
#
# figper = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, h_eff_change, '.', label="Effective Hamiltonian", color="red")
# plt.xlabel("time")
# plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
# plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()

############### PIMD for Fermions in different number of beads in a given Beta

# beads_array = np.array([3, 4, 5, 6, 8, 10])
print ("This is simulation is beta:3  itr: 4M  beads:  72")
beads_array = np.array([72])

e_tot_array = np.zeros(len(beads_array))
e_tot_stdv_array = np.zeros(len(beads_array))
ferm_energy_array = np.zeros(len(beads_array))
d_ferm_energy_array = np.zeros(len(beads_array))

for i, b in enumerate(beads_array):
    wp = math.sqrt(b) / (beta * hbar)
    start = time.time()
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change, s_o_dinom, s_o_num, hminush, kin_est1 = \
     langevin_dynamics_fermion_gaus(g_steps, dt, mass, beta, hbar, kboltz, w, wp, g, si, b, N_particles)
    stop = time.time()
    duration = stop - start
    print("Time of execution:", duration)
    print("Iteration:", i)

    number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=2000, data=temp_exp)
    number_of_blocks, kin_avg, stdv_kin = block_averaging(cutoff, block_size=4000, data=kin_est)
    number_of_blocks, potential_avg, stdv_potential = block_averaging(cutoff, block_size=4000, data=pot_est)
    number_of_blocks, avg_s_o_dinom, stdv_s_o_dinom = block_averaging(cutoff, block_size=4000, data=s_o_dinom)
    number_of_blocks, avg_s_o_num, stdv_s_o_num = block_averaging(cutoff, block_size=4000, data=s_o_num)
    number_of_blocks, avg_hminush, stdv_hminush = block_averaging(cutoff, block_size=4000, data=hminush)
    number_of_blocks, kin_avg1, stdv_kin1 = block_averaging(cutoff, block_size=4000, data=kin_est1)

    observable = (avg_s_o_num / avg_s_o_dinom)
    d_observable = (avg_s_o_num / avg_s_o_dinom) * np.sqrt((stdv_s_o_num / avg_s_o_num) ** 2 + (stdv_s_o_dinom / avg_s_o_dinom) ** 2)
    ferm_energy = observable + avg_hminush
    d_ferm_energy = math.sqrt(d_observable**2 + stdv_hminush**2)

# primitive vs virial
    print ("Kinetic Primitive: ", kin_avg1)
    print("Kinetic virial: ", kin_avg)

    print("Mean Temperature:", avg_temp_exp, "+-", stdv_temp)
    print("Fermion: Total Energy Estimator:", observable, "+-", d_observable)
    print("<H - H'>_H' :", avg_hminush, "+-", stdv_hminush)
    print("Fermionic Energy :", ferm_energy, "+-", d_ferm_energy)


    e_tot_array[i] = observable
    e_tot_stdv_array[i] = d_observable
    ferm_energy_array[i] = ferm_energy
    d_ferm_energy_array[i] = d_ferm_energy

print("Energy of H' :", e_tot_array, "+-", e_tot_stdv_array)
#BI Bogoliubov Inequality
print("Energy of H after BI:", ferm_energy_array, "+-", d_ferm_energy_array)

np.savez("2F_QHO_rep_gaus_beta_03_vs_beads_200M",
         ferm_energy_array=ferm_energy_array, d_ferm_energy_array=d_ferm_energy_array,
         e_tot_array=e_tot_array, e_tot_stdv_array=e_tot_stdv_array, beads_array=beads_array)

pass



figtotenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, e_tot_array, '.', label="Mean Total Energy", color="blue")
plt.errorbar(beads_array, e_tot_array, yerr=e_tot_stdv_array, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
plt.show()




 # Results for fermions

# beads_array = np.array([3, 4, 5, 6, 8, 10])
betas = np.array([1,  2,  3,  5, 10])
e_tot_f = np.array([])
stdv_f = np.array([])

e_tot_f = np.array([2.9524831472217326, 2.34320047, 2.24100858, 2.1173498253834837, 1.673302510683946])
stdv_f = np.array([0.10093912649307081, 0.10209355, 0.03397157, 0.030784098328776738, 0.22846661547608968])

figtotenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
q = np.linspace(0.2, 11, 1000)
p = ((hbar * w * (2 * np.exp(2 * q * hbar * w) + np.exp(q * hbar * w) + 1)) / (np.exp(2 * q * hbar * w) - 1))
plt.plot(q, p, 'g',  label="Analytical Result")
plt.plot(betas, e_tot_f, '.', label="PIMD", color="blue", linestyle='None')
plt.errorbar(betas, e_tot_f, yerr=stdv_f, ecolor="black", linestyle='None')
plt.xlabel("Beta")
plt.ylabel("Mean Total Energy")
plt.legend(loc='upper right')
plt.show()



# Fermions
# e05 =
# e1 =3.1488104538361865
# e2 = 2.5812701061054275
# e2_40M =2.423513606641101
# e3 =2.36505292171364
# e5 =2.3000630736334085
# e10 =2.032163420854567
# s05 =0.03364081299130081
# s1 =0.10038230645433016
# s2 = 0.10053454837907909
# s3 =0.04468156802744379
# s3_40M =0.03364081299130081
# s5 =0.030354862677618676
# s10 =0.22837578054151128
# E_b05=
# E_b1 = 2.9524831472217326
# E_b2 = 2.34320047
# E_b3 =2.18136903259555
# E_b3_40M =2.24100858
# E_b5 =2.1173498253834837
# E_10 =1.673302510683946
# d_E_b05=
# d_E_b1 = 0.10093912649307081
# d_E_b2 = 0.10209355
# d_E_b3 =0.0451619489281434
# d_E_b3_40M = 0.03397157
# d_E_b5 =0.030784098328776738
# d_E_10 =0.22846661547608968




# 2M 72 beads beta1
# Time of execution: 2255.669402360916
# Iteration: 0
# Mean Temperature: 0.9988394511702567 +- 0.0011123280949471204
# Fermion: Total Energy Estimator: 3.5134905674618775 +- 0.08896827158024327
# <H - H'>_H' : -0.19642332282404196 +- 0.01083358006727241
# Fermionic Energy : 3.3170672446378355 +- 0.08962544172861814



