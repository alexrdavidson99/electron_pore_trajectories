import numpy as np
import random
import matplotlib.pyplot as plt


my_nested_dict = {}

class Electrons(object):
    '''
    Using this class to represent electrons.
    '''

    def __init__(self, number, event_id, angle,start_position,end_position):
        self.number = number
        self.event_id = event_id
        self.angle = angle
        self.start_position = start_position
        self.end_position = end_position

    def __repr__(self):
        return (f"Electron(number={self.number}, event_id={self.event_id}, angle={self.angle}, "
                f"start_position={self.start_position}, end_position={self.end_position})")


def exponential_sequence(n, initial_value):
    sequence = [initial_value]
    total_e = []
    first_hit = False
    electron_list = []
    electron_counter = 0 # To keep track of unique electron numbers
    pt = 0
    start_position_counter = [0, 0, 0]

    for i in range(n - 1):

        if len(total_e) == 0 and not first_hit:
            energy = random.randint(150, 210)
            poisson_mean = 4.25 * (energy / 470) * np.exp(1 - (energy / 470))
            pos_0 = np.random.poisson(poisson_mean, 1)[0]
            next_number = pos_0
            my_nested_dict[i] = next_number
            first_hit = True
        else:
            next_number = sum(total_e)
            my_nested_dict[i] = next_number
            start_position_counter[0] += 1
            start_position_counter[1] += 1
            start_position_counter[2] += 1
            pt += 1
            print(f"test..{pt}")

        total_e.clear()
        for j in range(next_number):
            energy_2 = random.randint(150, 210)
            poisson_mean = 4.25 * (energy_2 / 470) * np.exp(1 - (energy_2 / 470))
            test = np.random.poisson(poisson_mean, 1)[0]

            # Create an Electron instance
            e_test = Electrons(electron_counter, i, 45, [0, 0, j], [start_position_counter[0], start_position_counter[1], start_position_counter[2]])
            electron_list.append(e_test)
            electron_counter += 1

            total_e.append(test)

        sequence.append(next_number)

    return sequence, electron_list

n = 6# Number of elements in the sequence
initial_value = 1

hist_g = []
all_electrons = []

for i in range(1):
    sequence, electrons = exponential_sequence(n, initial_value)
    if sequence is not None and electrons is not None:
        all_electrons.extend(electrons)
        y_hist = sequence[-1]
        hist_g.append(y_hist)

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

plt.show()
print(my_nested_dict)
# Print the list of all electrons for inspection
for electron in all_electrons:
    print(electron)
