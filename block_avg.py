
from MD_QHO_Functions import *

steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)


# block_size_array = np.array([5000])
block_size_array = np.linspace(1, 5000, 200).astype(int)
avg_array = np.zeros(len(block_size_array))
stdv_array = np.zeros(len(block_size_array))
number_of_blocks_array = np.zeros(len(block_size_array))

for i, value in enumerate(block_size_array):
    number_of_blocks, avg, stdv = block_averaging(cutoff, block_size=value, data=pot_est+kin_est)
    avg_array[i] = avg
    stdv_array[i] = stdv
    number_of_blocks_array[i] = number_of_blocks


figblock = plt.figure()
plt.plot(block_size_array, stdv_array, label="stdv", color="red")
plt.xlabel("Block Size")
plt.ylabel("STDV")
plt.legend()
plt.show()


np.savez("block_avg_pot_est", avg = avg_array, stdv = stdv_array, time=times, temp=temp_exp, pot=pot_est, kin=kin_est, heff=h_eff_change)
print("Average:", avg_array)
print ("stdv:", stdv_array)





