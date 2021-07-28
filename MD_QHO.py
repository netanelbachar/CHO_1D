from MD_QHO_Functions import *

# CONSTANTS
mass = 1
w = 1
hbar = 1
kboltz = 1
beta = 10

beads = 64

# Time
T_max = 2000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)


seed = 347885
np.random.seed(seed)

# Main
start = time.time()
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads)
stop = time.time()
duration = stop - start


print("Time of execution:", duration)
print("Mean Temperature:", np.mean(temp_exp[100000:]))
print("Mean Potential_estimator:", np.mean(pot_est[100000:]))
print("Mean Kinetic_estimator:", np.mean(kin_est[100000:]))


# Histograms
plt.figure(1)
n, bins, patches = plt.hist(pos, 50, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Position')
plt.ylabel('Counts')
plt.title('Position of Bead')
plt.grid(True)
# plt.show()

plt.figure(2)
n, bins, patches = plt.hist(vel, 50, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Velocity')
plt.ylabel('Counts')
plt.title('Velocity of Bead')
plt.grid(True)
# plt.show()


# Estimators

# figpotest = plt.figure()
# plt.plot(steps, pot_est, '.', label="P_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.plot(steps, kin_est, '.', label="K_Estimator vs step", color="blue")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()


# Conservation of Energy
print("Mean Total Energy:", np.mean(e_tot[100000:]))

figper = plt.figure()
plt.plot(times, e_change, '.', label="Energy % Change" , color="red")
plt.xlabel("time")
plt.ylabel("(Etot - E0)*100 / E0  [%]")
plt.ylim([-0.1, 0.1])
plt.legend()
# plt.show()
#
figenergy = plt.figure()
plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
plt.plot(times, potential, '.', label="Potential vs Time" , color="blue")
plt.xlabel("time")
plt.ylabel("Energy")
plt.legend()
# plt.show()

