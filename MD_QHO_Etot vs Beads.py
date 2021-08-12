from MD_QHO_Functions import *

# <Etot_est> vs Beads>

start = time.time()
mean_etot_array = np.zeros(len(beads_array))
stdv_etot_array = np.zeros(len(beads_array))

for i in range(0, len(beads_array)):
    seed = 347885
    np.random.seed(seed)
    print("beads", beads_array[i])
    steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
        langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array[i])
    tot_energy = kin_est + pot_est
    number_of_blocks, avg, stdv = block_averaging(cutoff, block_size=block_size, data=tot_energy)
    mean_etot_array[i] = avg
    stdv_etot_array[i] = stdv
stop = time.time()
duration = stop - start
print("Time of execution:", duration)

print("Mean ", mean_etot_array)
print("STDV ", stdv_etot_array)

np.savez("QHO_etotvsbeads_beta8_new_block", mean_e_tot_est=mean_etot_array, stdv=stdv_etot_array)

#QHO_etotvsbeads_beta10_new

figenest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, mean_etot_array, '.', label="Mean Total Energy Estimator", color="blue")
plt.errorbar(beads_array, mean_etot_array, yerr=stdv_etot_array, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
plt.show()


# Beta = 10   Bead =   32   Mean E = 0.49170702 +- 0.02319685
# Beta = 8    Bead =   32   Mean E = 0.49512283 +- 0.02536785
# Beta = 6    Bead =     Mean E =
# Beta = 4    Bead =      Mean E =
# Beta = 2    Bead =       Mean E =
# Beta = 1    Bead =       Mean E =

# Beta = 10   Bead = 20     Mean E = 0.4874987
# Beta = 8    Bead = 14     Mean E = 0.49260145
# Beta = 6    Bead = 14     Mean E = 0.48649026 +- 0.02031679
# Beta = 4    Bead = 10     Mean E = 0.49844144 +- 0.02061151    3M itr time execution 2683
# Beta = 2    Bead =  8     Mean E = 3M itr time execution 2683
# Beta = 1    Bead =  6     Mean E = 3M itr time execution 2683




