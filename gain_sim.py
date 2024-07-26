import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.constants

from energy_dis import accept_reject
from stepping import electron_yield, solve_for_intercept_time, step_position, step_energy
from consin_dis import cosine_dis


my_nested_dict = {}

class Electrons(object):
    '''
    Using this class to represent electrons.
    '''

    def __init__(self, number, event_id, energy,start_position,end_position):
        self.number = number
        self.event_id = event_id
        self.energy = energy
        self.start_position = start_position
        self.end_position = end_position

    def __repr__(self):
        return (f"Electron(number={self.number}, event_id={self.event_id}, energy={self.energy}, "
                f"start_position={self.start_position}, end_position={self.end_position})")


def exponential_sequence(n, initial_value):
    sequence = [initial_value]
    total_e = []
    first_hit = False
    electron_list = []
    yeild = []
    c = scipy.constants.speed_of_light * 1e-6 # Speed of light in mm/ns
    electron_counter = 0 # To keep track of unique electron numbers
    V = 1000
    m = 511e3
    d = 0.2
    l_over_d = 40
    r = (d / l_over_d) * 0.5
    E = V * (c ** 2) / (d * m)
    orientation = np.array([0., 0., 1])
    a0 = E * orientation # f=ma  where f = eE and a = E/m
    start_position = np.array([-r, 0, 0])
    initial_energy_distribution = accept_reject(1000) # Generate the energy values ejected from mcp

    for i in range(n - 1):

        if len(total_e) == 0 and not first_hit:
            energy = random.randint(250, 270)
            pos_0_yield = electron_yield(energy)
            yeild.append(pos_0_yield)
            next_number = pos_0_yield
            my_nested_dict[i] = next_number
            first_hit = True
        else:
            next_number = sum(total_e)
            my_nested_dict[i] = next_number

        total_e.clear()
        for j in range(next_number):
            initial_energy = random.choice(initial_energy_distribution)
            v_initial = np.sqrt((initial_energy / 511e3) * 2 * (c ** 2))
            vel_a, theta = cosine_dis(start_position, r)

            v1 = np.array([vel_a[0][0] * v_initial, vel_a[0][1] * v_initial, vel_a[0][2] * v_initial])
            ti = solve_for_intercept_time(start_position, v1, a0, r)
            print(f" time to hit {ti}")
            hit_position = step_position(start_position, v1, a0, ti)
            print(f"hit position {hit_position}")

            impact_energy = step_energy(v1, a0, ti, m)
            posssion_output = electron_yield(impact_energy)
            yeild.append(posssion_output)
            if hit_position[2] > d:
                print(f" start_position {start_position[2]/d} hit position is greater than the pore at location 3 {hit_position[2]/d}"
                      f" with hit energy {impact_energy}")
                break


            # Create an Electron instance
            e_test = Electrons(electron_counter, i, impact_energy, start_position.tolist(), hit_position.tolist())
            electron_list.append(e_test)
            electron_counter += 1
            start_position = hit_position
            total_e.append(posssion_output)



        sequence.append(next_number)

    return sequence, electron_list, yeild

n = 6# Number of elements in the sequence
initial_value = 1

hist_g = []
all_electrons = []
all_yields = []

for i in range(20):
    sequence, electrons, yeild = exponential_sequence(n, initial_value)
    if sequence is not None and electrons is not None:
        all_electrons.extend(electrons)
        y_hist = sequence[-1]
        hist_g.append(y_hist)
        all_yields.extend(yeild)

        y_plot = sequence
        x_plot = np.linspace(0, n, num=n)
        plt.plot(x_plot, y_plot, '-o')
    else:
        print(f"Iteration {i} returned None for sequence or electrons.")

# Plot the final results
plt.xlabel("Wall hits")
plt.ylabel("Number of electrons")

# Plot Poisson mean as a function of energy
plt.figure()
energy_test = np.linspace(1, 1000, 1000)
poisson_mean_test = 4.25 * (energy_test / 470) * np.exp(1 - (energy_test / 470))
plt.plot(energy_test, poisson_mean_test)

plt.xlabel("Energy")
plt.ylabel("Poisson Mean")
plt.figure()
energies = [electron.energy for electron in all_electrons]
plt.hist(energies, bins=100, alpha=0.6, label='energy histogram')

plt.show()
print(my_nested_dict)
# Print the list of all electrons for inspection
#for electron in all_electrons:
   # print(electron)
#print(np.mean(all_yields))
