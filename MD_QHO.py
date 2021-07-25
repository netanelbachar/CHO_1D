from MD_QHO_Functions import *

# CONSTANTS
mass = 1
w = 1
hbar = 1
kboltz = 1
beta = 10

beads = 16

# Time
T_max = 5000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)


seed = 347885
np.random.seed(seed)

# Main
start = time.time()
steps, times, kin, potential, e_tot, temp_exp, pot_est, kin_est = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)
stop = time.time()
duration = stop - start


print("Time of execution:", duration)
print("Mean Temperature:", np.mean(temp_exp[100000:]))
print("Mean Potential_estimator:", np.mean(pot_est[250000:]))
print("Mean Kinetic_estimator:", np.mean(kin_est[250000:]))


# figpotest = plt.figure()
# plt.plot(steps, pot_est, '.', label="P_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()

# figpotest = plt.figure()
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()

# fig2 = plt.figure()
# plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
# plt.plot(times, potential, '.', label="Potential vs Time" , color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()