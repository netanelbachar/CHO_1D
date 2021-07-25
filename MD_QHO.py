from MD_QHO_Functions import *

# CONSTANTS
mass = 1
w = 1
hbar = 1
kboltz = 1
beta = 10

beads = 32

# Time
T_max = 2000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)


seed = 347885
np.random.seed(seed)

# Main
start = time.time()
steps, times, kin, e_tot, temp_exp, pot_est, kin_est = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)
stop = time.time()
duration = stop - start


print("Time of execution:", duration)
print("Mean Temperature:", np.mean(temp_exp[100000:]))
print("Mean Potential_estimator:", np.mean(pot_est[100000:]))
print("Mean Kinetic_estimator:", np.mean(kin_est[100000:]))


# figpotest = plt.figure()
# plt.plot(steps, pot_est, '.', label="P_Estimator vs time", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()

