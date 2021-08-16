
from MD_QHO_Functions import *


# Calculated Mean Estimators (kinetic or potential) for a given number of beads

start = time.time()
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)
number_of_blocks, avg, stdv = block_averaging(cutoff, block_size=block_size, data=(kin_est+pot_est))
stop = time.time()
duration = stop - start

print("Time of execution:", duration)
print("Mean Temperature:", np.mean(temp_exp[cutoff:]))
print("Mean Potential_estimator:", np.mean(pot_est[cutoff:]))
print("Mean Kinetic_estimator:", np.mean(kin_est[cutoff:]))
print("Mean Total Energy Estimator:", np.mean(kin_est[cutoff:]) + np.mean(pot_est[cutoff:]))
print("BA_Mean Total Energy Estimator:", avg, "+-", stdv)


# Histograms

# plt.figure(1)
# plt.rcParams.update({'font.size': 13})
# n, bins, patches = plt.hist(pos, 50, density=True, facecolor='b', alpha=0.75)
# plt.xlabel('Position')
# plt.ylabel('Counts')
# plt.title('Position of Bead')
# plt.xlim([-1, 1])
# plt.grid(True)
# plt.show()

plt.figure(2)
plt.rcParams.update({'font.size': 13})
n, bins, patches = plt.hist(vel, 50, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Velocity')
plt.ylabel('Counts')
plt.title('Velocity of Bead')
plt.xlim([-1.2, 1.2])
plt.grid(True)
plt.show()
#
# # Conservation of Energy and Temperature
#
# print("Mean Total Energy:", np.mean(e_tot[cutoff:]))
#
# figtemp = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, temp_exp, '.', label="Temperature", color="red")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()
#
figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, e_change, '.', label="Energy % Change", color="red")
plt.xlabel("time")
plt.ylabel("(E(t) - E0)*100 / E0  [%]")
plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()
#
# figenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
# plt.plot(times, potential, '.', label="Potential vs Time" , color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()
#
# # Estimators
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
#
#
# # Effective Hemiltonian for Energy Conservation
#
# figper = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, h_eff_change, '.', label="Effective Hemiltonian", color="red")
# plt.xlabel("time")
# plt.ylabel("(H_eff(t) - H_eff(0))*100 / H_eff(0) [%]")
# plt.ylim([-0.1, 0.1])
# plt.legend()
# plt.show()

np.savez("QHO_beta10", avg=avg, stdv=stdv, time=times, pos=pos, vel=vel, kin=kin, potential= potential, e_tot=e_tot,
         e_change=e_change, temp_exp=temp_exp, pot_est=pot_est, kin_est=kin_est, h_eff_change=h_eff_change)

