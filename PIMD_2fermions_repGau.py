from MD_QHO_Functions_2bosons import *

start = time.time()
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est, kin_est, h_eff_change = \
 langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, wp, g, si, beads, N_particles)
stop = time.time()
duration = stop - start
print("Time of execution:", duration)

print("bla", len(e_tot))
figenergy = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
plt.plot(times, potential, '.', label="Potential vs Time", color="blue")
plt.xlabel("time")
plt.ylabel("Energy")
plt.legend()
plt.show()


figper = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, e_change, '.', label="% E_tot change", color="red")
plt.xlabel("time")
plt.ylabel("(E(t) - E(0))*100 / E(0) [%]")
plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()

# figtemp = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, temp_exp, '.', label="Temperature", color="red")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()

plt.figure(1)
plt.rcParams.update({'font.size': 13})
n, bins, patches = plt.hist(pos, 100, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Position')
plt.ylabel('Counts')
plt.title('Position of Bead')
# plt.xlim([-1, 1])
plt.grid(True)
plt.show()

plt.figure(2)
plt.rcParams.update({'font.size': 13})
n, bins, patches = plt.hist(vel, 100, density=True, facecolor='b', alpha=0.75)
plt.xlabel('Velocity')
plt.ylabel('Counts')
plt.title('Velocity of Bead')
# plt.xlim([-1.2, 1.2])
plt.grid(True)
plt.show()

