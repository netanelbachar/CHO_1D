from MD_QHO_Functions import *

# What is left: Run this code for beta = 6,4,2,1 in dt = 10-**-3

# <Etot_est> vs Beads>

start = time.time()
mean_e_tot_est, stdv = \
    langevin_dynamics_beads(g_steps, cutoff, num_blocks, dt, mass, beta, hbar, kboltz, w, beads_array)
stop = time.time()
duration = stop - start

print("Time of execution:", duration)
print("Mean ", mean_e_tot_est)
print("STDV", stdv)

figenest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(beads_array, mean_e_tot_est, '.', label="Mean Total Energy Estimator", color="blue")
plt.errorbar(beads_array, mean_e_tot_est, yerr=stdv, ecolor="black")
plt.xlabel("Beads")
plt.ylabel("Mean Total Energy")
plt.legend(loc='lower right')
plt.show()

# Beta = 10   Bead = 20     Mean E = 0.4874987
# Beta = 8    Bead = 14     Mean E = 0.49260145
# Beta = 6    Bead = 14     Mean E = 0.48649026 +- 0.02031679
# Beta = 4    Bead = 10     Mean E = 0.49844144 +- 0.02061151    3M itr time execution 2683
# Beta = 2    Bead =  8     Mean E = 3M itr time execution 2683
# Beta = 1    Bead =  6     Mean E = 3M itr time execution 2683




