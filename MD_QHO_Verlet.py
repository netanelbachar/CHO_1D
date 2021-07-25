# Quantum Harmonic Oscillator in One Dimension

# IMPORTS
from functions_QHO import *

# CONSTANTS
N_particles = 1
mass = 1
w = 1
kboltz = 1
hbar = 1
temperature = 1
beta = 1 / (temperature * kboltz)
beads = 1
wp = math.sqrt(beads) / (beta * hbar)
# Time

T_max = 20
dt = 1 * 10 ** (-4)
g_steps = int(T_max / dt)
print("Number of Steps:", g_steps)
seed = 347885
np.random.seed(seed)


# Main Code

steps, times, pos, vx_n, kin, pot, e_tot, percent_change, temp_exp = verlet_propagator(g_steps, dt, mass, w, beta, hbar, kboltz, wp, beads)


num = int(g_steps * 0.1)
print("Mean Energy:", sum(e_tot[num:])/len(e_tot[num:]))


# PLOTS

fig2 = plt.figure()
plt.plot(times, e_tot, '.', label="Total Energy vs Time" , color="black")
plt.plot(times, kin, '.', label="Kinetic vs Time" , color="red")
plt.plot(times, pot, '.', label="Potential vs Time" , color="blue")
plt.xlabel("time")
plt.ylabel("Energy")
plt.legend()
plt.show()

figper = plt.figure()
plt.plot(times, percent_change, label="Percent Change", color="red")
plt.xlabel("dt")
plt.ylabel("% Change ((Ei-E0)*100/E0) ")
plt.ylim([-0.1, 0.1])
plt.legend()
plt.show()


fig1 = plt.figure()
plt.plot(times, pos, '.', label="Position of Particle" , color="black")
plt.xlabel("time")
plt.legend()
plt.show()

plt.figure(figsize=[10, 8])
n, bins, patches = plt.hist(x=pos, bins=30, density=True, alpha=0.7,rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Position', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.title('MD_CHO NVE : Position Histogram', fontsize=15)
plt.show()

# figtemp = plt.figure()
# plt.plot(times, temp_exp, '.', label="Temperature vs time" , color="black")
# plt.xlabel("time")
# plt.ylabel("Temperature")
# plt.legend()
# plt.show()



