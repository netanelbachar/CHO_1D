from MD_QHO_Functions_2bosons import *


# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change, s_o_dinom, s_o_num = \
#     langevin_dynamics_fermion(g_steps, dt, mass, beta, hbar, kboltz, w, wp, beads, N_particles)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=2000, data=temp_exp)
# number_of_blocks, kin_avg, stdv_kin = block_averaging(cutoff, block_size=1000, data=kin_est)
# number_of_blocks, potential_avg, stdv_potential = block_averaging(cutoff, block_size=1000, data=pot_est)
# number_of_blocks, s_o_dinom, stdv_s_o_dinom = block_averaging(cutoff, block_size=1000, data=s_o_dinom)
# number_of_blocks, s_o_num, stdv_s_o_num = block_averaging(cutoff, block_size=1000, data=s_o_num)
#
# observable = (s_o_num / s_o_dinom)
# d_observable = (s_o_num / s_o_dinom) * np.sqrt((stdv_s_o_num / s_o_num)**2 + (stdv_s_o_dinom / s_o_dinom)**2)
#
# print("Mean classical Temperature:", avg_temp_exp, "+-", stdv_temp)
# print("Potential Estimator-Fermion:", observable, "+-", d_observable)
# print("Total Energy Estimator-Fermion:", 2 * observable, "+-", 2 * d_observable)
#
#
# np.savez("2F_QHO_beta10", observable=2*observable, d_observable=2*d_observable)


beads_array = np.array([3, 4, 5, 6, 8, 10])
e_tot_array = np.zeros(len(beads_array))
e_tot_stdv_array = np.zeros(len(beads_array))

for i, b in enumerate(beads_array):
    wp = math.sqrt(b) / (beta * hbar)

    start = time.time()
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change, s_o_dinom, s_o_num = \
        langevin_dynamics_fermion(g_steps, dt, mass, beta, hbar, kboltz, w, wp, b, N_particles)
    stop = time.time()
    duration = stop - start
    print("Time of execution:", duration)
    print("Iteration:", i)

    number_of_blocks, avg_temp_exp, stdv_temp = block_averaging(cutoff, block_size=2000, data=temp_exp)
    number_of_blocks, kin_avg, stdv_kin = block_averaging(cutoff, block_size=1000, data=kin_est)
    number_of_blocks, potential_avg, stdv_potential = block_averaging(cutoff, block_size=1000, data=pot_est)
    number_of_blocks, avg_s_o_dinom, stdv_s_o_dinom = block_averaging(cutoff, block_size=4000, data=s_o_dinom)
    number_of_blocks, avg_s_o_num, stdv_s_o_num = block_averaging(cutoff, block_size=4000, data=s_o_num)

    observable = (avg_s_o_num / avg_s_o_dinom)
    d_observable = (avg_s_o_num / avg_s_o_dinom) * np.sqrt((stdv_s_o_num / avg_s_o_num) ** 2 + (stdv_s_o_dinom / avg_s_o_dinom) ** 2)

    print("Mean classical Temperature:", avg_temp_exp, "+-", stdv_temp)
    print("Potential Estimator-Fermion:", observable, "+-", d_observable)
    print("Total Energy Estimator-Fermion:", 2 * observable, "+-", 2 * d_observable) # times 2 since kin = pot

    e_tot_array[i] = 2 * observable
    e_tot_stdv_array[i] = 2 * d_observable

print("e array:", e_tot_array)
print("e stdv array:", e_tot_stdv_array)

np.savez("2F_QHO_beta03_vs_beads_400K",  e_tot_array=e_tot_array, e_tot_stdv_array=e_tot_stdv_array, beads_array=beads_array)

pass


figtotenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, e_tot_array, '.', label="Mean Total Energy", color="blue")
plt.errorbar(beads_array, e_tot_array, yerr=e_tot_stdv_array, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
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

figtemp = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, temp_exp, '.', label="Temperature", color="red")
plt.xlabel("time")
plt.ylabel("Temperature")
plt.legend()
plt.show()

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
# figkinest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()



 # Results for fermions


betas = np.array([0.3, 0.5, 0.8, 1,  2,  3,  5, 10])
e_tot_f = np.array([7.3768407, 4.758762338333333, 3.350115618333333, 2.8454311433333337, 2.189582618333333, 2.0275577000000005, 1.5472921679999998, 1.1172278766666668])
stdv_f = np.array([0.16043649166666665, 0.12524463500000002, 0.11092883333333332, 0.16780319000000002, 0.26941090333333334, 0.6055924033333333, 1.557349254, 2.0800425866666665])

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
e03_400 = [7.1507412,  7.48628518, 7.41496091, 7.48702217, 7.34209828, 7.37993646]
e05_200 = [4.5357815, 4.71327161, 4.63188427, 4.78677412, 4.56457624, 4.5373544]
e_08_400 = [3.25806633, 3.34994199, 3.36081533, 3.38470269, 3.42320071, 3.32396666]
e05_400 = [4.60585359, 4.82351404, 4.77437378, 4.79854493, 4.75275367, 4.79753402]
e1 = [2.84610148, 2.89407776, 2.92581561, 2.89597672, 2.75570742, 2.75490787]
e2 = [2.10056319, 2.27104433, 2.15804499, 2.18094994, 2.20390157, 2.22299169]
e2_400 = [2.08649002, 2.21690686, 2.14941599, 2.28565752, 2.2159583,  2.11484461]
e3 = [1.86564714, 2.11722812, 1.88244831, 1.81771985, 2.22110619, 2.26119659]
e5 = [1.68192383, 1.79057334, 1.64589027,  1.34222339,  1.27585001, -5.12486395]
e10 = [0.98475069, 1.1838517,  1.26417233, 0.82768363, 1.4252744,  1.01763451]
s03_400 = [0.12828373, 0.14530887, 0.15358778, 0.16500423, 0.17702083, 0.19341351]
s05_200 = [0.13980038, 0.15962774, 0.17167668, 0.1948065,  0.19701865, 0.22425398]
s05_400 = [0.10151122, 0.10789534, 0.11992312, 0.12702102, 0.13992461, 0.1551925 ]
s_08_400 = [0.0869377,  0.09853181, 0.1054197,  0.11335768, 0.12757182, 0.13375429]
s1 = [0.12451232, 0.13613391, 0.17023989, 0.16922247, 0.18744647, 0.21926408]
s2 = [0.18170543, 0.20825258, 0.23467584, 0.23030649, 0.34495139, 0.41657369]
s2_400 = [0.12191641, 0.14630187, 0.1636521,  0.18638925, 0.20479888, 0.19709186]
s3 = [0.32504659, 0.35975937, 0.39738406, 0.33755085, 1.07726411, 1.13654944]
s5 = [1.76264082,   0.75965325,   1.50750499,   0.6326556,    3.12429161, -53.92599749]
s10 = [1.0709464,  0.8940214,  5.84186777, 0.86974325, 1.91786077, 1.88581593]


