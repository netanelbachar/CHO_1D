from MD_QHO_Functions_2bosons import *


def oscillator_potential(mass, w, x):
    potential = (0.5 * mass * w ** 2 * x ** 2).sum()
    return potential


def springs_potential(mass, wp, x):
    potential = 0
    for j in range(len(x)):
        if j == (len(x) - 1):
            potential += 0.5 * mass * wp ** 2 * (x[0] - x[j]) ** 2
        else:
            potential += 0.5 * mass * wp ** 2 * (x[j + 1] - x[j]) ** 2
    return potential


def oscillator_force(mass, w, beads, x):
    f1 = - mass * w ** 2 * x / beads
    return f1


def ring_springs_force(mass, x):
    wp = math.sqrt(len(x)) / (beta * hbar)
    c2 = - mass * wp ** 2
    f2 = np.zeros(len(x))
    for j in range(len(x)):
        if j == (len(x) - 1):
            f2[j] = c2 * (2 * x[j] - x[0] - x[j - 1])
        elif j == 0:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[-1])
        else:
            f2[j] = c2 * (2 * x[j] - x[j + 1] - x[j - 1])
    return f2



beads_n = beads // N_particles
b = 4
N = 2

x = initial_position(mass, w, hbar, (N, b))
vx = initial_velocity(mass, beta, beads)
x_long = np.concatenate((x[0], x[1]),  axis=None)

f0 = ring_springs_force(mass, x[0]) + oscillator_force(mass, w, b, x[0])
f1 = ring_springs_force(mass, x[1]) + oscillator_force(mass, w, b, x[1])
f12 = ring_springs_force(mass, x_long) + oscillator_force(mass, w, b, x_long) #  not good

print ("f0", f0)
print ("f1", f1)
print ("f12", f12)































# v_oo = springs_potential(mass, wp, x1) + springs_potential(mass, wp, x2) + \
#        (oscillator_potential(mass, w, x1) + oscillator_potential(mass, w, x2)) / beads_n
# v_o = springs_potential(mass, wp, x) + oscillator_potential(mass, w, x)
# exp_v_oo = np.exp(- beta * v_oo)
# exp_v_o = np.exp(- beta * v_o)
