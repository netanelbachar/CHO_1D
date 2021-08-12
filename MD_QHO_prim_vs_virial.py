from MD_QHO_Functions import *

# See the Fluctuations in the Potential Estimator and Kinetic Estimator  vs time with different number of Beads


# Primitive

# start = time.time()
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est0, kin_est0, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[0])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est1, kin_est1, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[1])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est2, kin_est2, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[2])
# steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est3, kin_est3, h_eff_change = \
#     langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[3])
# stop = time.time()
# duration = stop - start
#
# np.savez("QHO_estimator_primitive", pot_est0=pot_est0, kin_est0=kin_est0, pot_est1=pot_est1, kin_est1=kin_est1
#          , pot_est2=pot_est2, kin_est2=kin_est2, pot_est3=pot_est3, kin_est3=kin_est3)
#
# figkinest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, kin_est0, '.', label="Kin_Est 10 Beads", color="black")
# plt.plot(times, kin_est1, '.', label="Kin_Est 20 Beads", color="red")
# plt.plot(times, kin_est2, '.', label="Kin_Est 30 Beads", color="blue")
# plt.plot(times, kin_est3, '.', label="Kin_Est 40 Beads", color="green")
# plt.xlabel("time")
# plt.ylabel("Kinetic Energy Estimator")
# plt.legend()
# plt.show()
#
# figpotest = plt.figure()
# plt.rcParams.update({'font.size': 13})
# plt.plot(times, pot_est0, '.', label="Pot_Est 10 Beads", color="black")
# plt.plot(times, pot_est1, '.', label="Pot_Est 20 Beads", color="red")
# plt.plot(times, pot_est2, '.', label="Pot_Est 30 Beads", color="blue")
# plt.plot(times, pot_est3, '.', label="Pot_Est 40 Beads", color="green")
# plt.xlabel("time")
# plt.ylabel("Potential Energy Estimator")
# plt.legend()
# plt.show()

# Virial

start = time.time()
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est0, kin_est0, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[0])
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est1, kin_est1, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[1])
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est2, kin_est2, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[2])
steps, times, pos, vel, kin, potential, e_tot, e_change, temp_exp, pot_est3, kin_est3, h_eff_change = \
    langevin_dynamics(g_steps, dt, mass, beta, hbar, kboltz, w, beads_array_p[3])
stop = time.time()
duration = stop - start

np.savez("QHO_estimator_virial", pot_est0=pot_est0, kin_est0=kin_est0, pot_est1=pot_est1, kin_est1=kin_est1
         , pot_est2=pot_est2, kin_est2=kin_est2, pot_est3=pot_est3, kin_est3=kin_est3)

figkinest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, kin_est0, '.', label="Kin_Est 10 Beads", color="black")
plt.plot(times, kin_est1, '.', label="Kin_Est 20 Beads", color="red")
plt.plot(times, kin_est2, '.', label="Kin_Est 30 Beads", color="blue")
plt.plot(times, kin_est3, '.', label="Kin_Est 40 Beads", color="green")
plt.xlabel("time")
plt.ylabel("Kinetic Energy Estimator Virial")
plt.legend()
plt.show()

figpotest = plt.figure()
plt.rcParams.update({'font.size': 13})
plt.plot(times, pot_est0, '.', label="Pot_Est 10 Beads", color="black")
plt.plot(times, pot_est1, '.', label="Pot_Est 20 Beads", color="red")
plt.plot(times, pot_est2, '.', label="Pot_Est 30 Beads", color="blue")
plt.plot(times, pot_est3, '.', label="Pot_Est 40 Beads", color="green")
plt.xlabel("time")
plt.ylabel("Potential Energy Estimator")
plt.legend()
plt.show()




