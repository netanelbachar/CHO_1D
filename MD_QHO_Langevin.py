# Quantum Harmonic Oscillator in One Dimension for one particle
# I am trying to get the estimator vs time to converge

# IMPORTS


from functions_QHO import *

# CONSTANTS
seed = 347885
np.random.seed(seed)
mass = 1
w = 1
kboltz = 1
hbar = 1
temperature = 0.1
beta_s = 1 / (temperature * kboltz)
# temperature = np.array([10, 5, 1, 0.2, 0.1])  # 298.15  # [K]

beads = 8

# Time
T_max = 2000
dt = 1 * 10 ** (-2)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)


#Main Code 1

wp = math.sqrt(beads) / (beta_s * hbar)
start = time.time()
steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp, pot_est, kin_est = \
    langevin_thermostat(g_steps, dt, mass, w, beta_s, hbar, kboltz, wp, beads)
stop = time.time()
duration = stop - start

tot_est = kin_est + pot_est

print("Time of execution:", duration)
print ("Mean Temperature:", np.mean(temp_exp))
print ("Mean Potential_estimator:", np.mean(pot_est[100000:]))

plt.figure(figsize=[10, 8])
n, bins, patches = plt.hist(x=pos, bins=30, density=True, alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Position', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.title('MD_CHO NVT : Position Histogram', fontsize=15)
plt.show()

plt.figure(figsize=[10, 8])
n, bins, patches = plt.hist(x=vx_n, bins=30, density=True, alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Velocity', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.title('MD_CHO NVT : Velocity Histogram', fontsize=15)
plt.show()


fig1 = plt.figure()
plt.plot(times, pos, '.', label="Position of Particle" , color="black")
plt.xlabel("time")
plt.legend()
plt.show()

figtemp = plt.figure()
plt.plot(times, temp_exp, '.', label="Temperature vs time" , color="black")
plt.xlabel("time")
plt.ylabel("Temperature")
plt.legend()
plt.show()

figpotest = plt.figure()
plt.plot(steps, pot_est, '.', label="P_Estimator vs time", color="blue")
plt.xlabel("Steps")
plt.ylabel("Potential Estimator")
plt.legend()
plt.show()

#Main Code 2

# start = time.time()
# beads_array = np.array([1, 2, 4, 8, 16, 24, 32, 64])
# pot_est_array_s = estimator_potential_array(g_steps, dt, mass, w, beta_s, hbar, kboltz, beads_array)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# figpotest = plt.figure()
# plt.plot(beads_array, pot_est_array_s, label="Pot_Estimator vs Beads", color="black")
# plt.xlabel("P")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()

####
#Plot
#




# figkinest = plt.figure()
# plt.plot(steps, kin_est, '.', label="K_Estimator vs time", color="red")
# plt.xlabel("Steps")
# plt.ylabel("Kinetic Estimator")
# plt.legend()
# plt.show()
#
# figtotest = plt.figure()
# plt.plot(steps, tot_est, '.', label="E_Estimator vs time", color="black")
# plt.xlabel("Steps")
# plt.ylabel("Total Energy Estimator")
# plt.ylim([0, 11])
# plt.legend()
# plt.show()
####

# plt.figure(figsize=[10, 8])
# n, bins, patches = plt.hist(x=pos, bins=50, density=True, alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('Position',fontsize=15)
# plt.ylabel('Count',fontsize=15)
# plt.title('MD_CHO NVT : Position Histogram',fontsize=15)
# plt.show()
#
#
# plt.figure(figsize=[10, 8])
# n, bins, patches = plt.hist(x=vx_n, bins=50, density=True, alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('Velocity', fontsize=15)
# plt.ylabel('Count', fontsize=15)
# plt.title('MD_CHO NVT : Velocity Histogram', fontsize=15)
# plt.show()

##### When estimator works then
# beta1 = 1 / (temperature[0] * kboltz)
# beta2 = 1 / (temperature[1] * kboltz)
# beta3 = 1 / (temperature[2] * kboltz)
# beta4 = 1 / (temperature[3] * kboltz)
# beta5 = 1 / (temperature[4] * kboltz)
# beta = np.array([beta1, beta2, beta3, beta4, beta5])

# start = time.time()
# beads_array = np.array([1, 2, 4, 8, 16, 24, 32, 64, 72])
# pot_est_array_s = estimator_potential_array(g_steps, dt, mass, w, beta_s, hbar, kboltz, beads_array)
# stop = time.time()
# duration = stop - start
# print("Time of execution:", duration)
#
# figpotest = plt.figure()
# plt.plot(beads_array, pot_est_array_s, label="Pot_Estimator vs Beads", color="black")
# plt.xlabel("P")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()

# pot_est_array1 = estimator_potential_array(g_steps, dt, mass, w, beta1, hbar, kboltz, beads_array)
# pot_est_array2 = estimator_potential_array(g_steps, dt, mass, w, beta2, hbar, kboltz, beads_array)
# pot_est_array3 = estimator_potential_array(g_steps, dt, mass, w, beta3, hbar, kboltz, beads_array)
# pot_est_array4 = estimator_potential_array(g_steps, dt, mass, w, beta4, hbar, kboltz, beads_array)
# pot_est_array5 = estimator_potential_array(g_steps, dt, mass, w, beta5, hbar, kboltz, beads_array)
#
# pot_est_array = np.array([pot_est_array1[-1], pot_est_array2[-1], pot_est_array3[-1], pot_est_array4[-1], pot_est_array5[-1]])
#
# # PLOTS NVT
#
# figpotest = plt.figure()
# plt.plot(beads_array, pot_est_array1, label="Beta = 0.1 ", color="black")
# plt.plot(beads_array, pot_est_array2, label="Beta = 0.5", color="blue")
# plt.plot(beads_array, pot_est_array3, label="Beta = 1", color="red")
# plt.plot(beads_array, pot_est_array4, label="Beta = 5", color="green")
# plt.plot(beads_array, pot_est_array5, label="Beta = 10", color="yellow")
# plt.xlabel("number of beads")
# plt.ylabel("Potential Estimator")
# plt.ylim([0, 30])
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.plot(beta, pot_est_array, label="Pot_Estimator vs Beta", color="black")
# plt.xlabel("Beta")
# plt.ylabel("Potential Estimator")
# plt.legend()
# plt.show()




###########

# fig2 = plt.figure()
# plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
# plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
# plt.plot(times, pot, '.', label="Potential vs Time" , color="blue")
# plt.xlabel("time")
# plt.ylabel("Energy")
# plt.legend()
# plt.show()

# fig1 = plt.figure()
# plt.plot(times, temp_exp, '.', label="Temperature vs time" , color="black")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()
#
# figper = plt.figure()
# plt.plot(times, percent_change, label="Percent Change", color="red")
# plt.xlabel("dt")
# plt.ylabel("% Change ((Ei-E0)*100/E0) ")
# plt.ylim([-0.5, 0.5])
# plt.legend()
# plt.show()

