import numpy
import numpy as np

from MD_QHO_Functions import *

# What is left: Check that it works 100%

num_of_blocks_array = np.array([3, 6, 8, 10, 14, 18, 26, 30, 32, 38, 46, 56, 66, 88, 102, 124, 148, 200, 250])
# num_of_blocks_array = np.array([3, 6, 8, 10, 14, 18, 26, 30, 32, 38])
avg_array = np.zeros(len(num_of_blocks_array))
stdv_array = np.zeros(len(num_of_blocks_array))

steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)

for i, value in enumerate(num_of_blocks_array):
    avg, stdv = block_averaging(cutoff, num_of_blocks=value, data=kin_est)
    avg_array[i] = avg
    stdv_array[i] = stdv

figblock = plt.figure()
plt.plot(num_of_blocks_array, stdv_array, label="stdv", color="red")
plt.xlabel("Number of Blocks")
plt.ylabel("STDV")
plt.legend()
plt.show()

numpy.savez("block_avg_file", temp=temp_exp, pot=pot_est, kin=kin_est, heff=h_eff_change)

print("Average:", avg_array)

print ("stdv:", stdv_array)
