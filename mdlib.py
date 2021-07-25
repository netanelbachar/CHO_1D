# This is a Library created by Netanel Bachar Schwartz of all usefull functions for MD simulations



def initial_position_1d():
    sp.random.seed(357922)   # To obtain the same random numbers (each number in () creates different cenario
    x = np.random.uniform(0, 1)
    return x


def initial_velocity_1d(number_of_particles, temp):
    sigma = np.sqrt(kboltz * temp / mass)  # -x^2/2s^2 = -v^2/(2kT/m) -> s^2 = kT/m
    vx = np.random.normal(0.0, sigma, number_of_particles)
    # Velocity of Center of Mass  = 0
    vx = vx - (mass * vx)/(mass * number_of_particles)      # only because all masses are equal!
    return vx, sigma


def harmonic_potential_1d(mass, w, x):
    potential = 0.5 * mass * w**2 * x**2
    return potential


def kinetic_energy_1d(mass, vx):
    kinetic_x = 0.5 * mass * vx ** 2
    return kinetic_x


def langevin_gamma(dt):  # This is the friction term
    gamma = 0.01 / dt
    return gamma


def langevin_xsi():  # This is the noise term
    xsi = np.random.normal(0, 1, 1)
    return xsi


def langevin_dynamics_1d(v, time_step):
    vel = math.exp(-1 * l_gamma(time_step) * time_step / 2) * v + \
                (1 / (math.sqrt(mass * beta))) * math.sqrt(1-math.exp(-1 * l_gamma(time_step) * time_step)) * l_xsi()
    return vel

