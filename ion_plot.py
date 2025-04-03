import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import pandas as pd



def solve_for_intercept_time(x0, v0, acc, target_distance):
    """
    Solve for the intercept time when the particle reaches the target distance in the y direction.
    """
    # Define the polynomial coefficients for the equation of motion in the y direction
    coeffs = [
        0.5 * acc[1],  # t^2 term
        v0[1],  # t term
        x0[1] - target_distance  # constant term
    ]

    # Solve the quadratic equation
    roots = np.roots(coeffs)

    # Filter out complex roots and negative times
    real_roots = roots[np.isreal(roots) & (roots >= 0)]

    if len(real_roots) == 0:
        raise ValueError("No valid intercept time found.")

    # Return the smallest positive real root
    return np.min(real_roots)


def step_position(x0, v0, acc, time):
    """
    Position after time step
    """
    x0 = np.array(x0)
    v0 = np.array(v0)
    acc = np.array(acc)
    return x0 + v0 * time + 0.5 * acc * time ** 2


def rotate_coordinates(x, y, theta):
    """
    Rotate coordinates
    """
    angle_degrees = - theta
    # Convert the angle to radians
    angle_radians = np.deg2rad(angle_degrees)

    # Create the rotation matrix
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians), np.cos(angle_radians)]
    ])
    coordinates = np.stack((x, y))
    # Apply the rotation matrix to the coordinates
    rotated_coordinates = np.dot(rotation_matrix, coordinates)
    x_rotate, y_rotate = rotated_coordinates

    return x_rotate, y_rotate


def electron_trajectory(target_distance, q_t, mass, initial_energy, theta_deg, voltage):
    """
    Simulates the trajectory of an electron (or ion) in an electric field until it reaches the target distance in
    the z-axis.
    
    Parameters:
    - initial_energy: Initial kinetic energy in Joules
    - target_distance: Target distance in meters for the z-axis
    - mass: Mass of the particle in kg
    - theta_deg: Electric field angle in degrees (default 8)
    - dt: Time step in seconds (default 1e-12)
    
    Returns:
    - t: Time array
    - x, y, z: Position arrays
    """
    c = scipy.constants.speed_of_light
    print(f"mass {mass} eV")
    v_i = np.sqrt((2*initial_energy)/mass)*c
    print(f"initial velocity {v_i/c} m/s")

    theta = np.deg2rad(theta_deg)
    cathode_field = voltage / target_distance
    # Calculate acceleration components
    a_y = cathode_field * q_t
    print(f"acceleration in y direction {a_y} m/s^2")

    # Initialize lists for time, position, and velocity

    x = []
    y = []

    vx = [v_i*np.sin(theta)]
    vy = [v_i * np.cos(theta)]

    # Solve for intercept time
    intercept_time = solve_for_intercept_time([0, 0], [vx[0], vy[0]], [0, a_y], target_distance)
    t = [intercept_time]
    print(f"Time to intercept target {intercept_time*1e9} ns")
    x_end = step_position([0, 0], [vx[0], vy[0]], [0, a_y], intercept_time)
    print(f"end position {x_end}")
    dt = 1000
    # Simulate the trajectory
    for i in range(dt):
        x_end = step_position([0, 0], [vx[0], vy[0]], [0, a_y], intercept_time*(i/dt))
        x.append(x_end[0])
        y.append(x_end[1])
    return x, y, t

def from_csv_energy_value_cal(pore_length, sample_start_position):
    random_value = np.random.choice(sample_start_position["start_position"])
    #cal_random_value = pore_length - random_value*1e-3
    cal_random_value = random_value*1e-3
    print(f"random {random_value}")
    return 750*cal_random_value/pore_length

# Parameters
target_distance = 1550e-6  # Target distance in meters
pore_length = 460e-6  # Pore length in meters   
Q_t = 1.76e11  # Mass of electron in kg
intital_energy_ion = 0
#initial_energies = np.linspace(intital_energy_ion, 700, 500)

sample_start_position = pd.read_csv("sampled_start_positions_ar.csv", header=None, names=["start_position"], skiprows=1)
sampled_energy = [0] + [from_csv_energy_value_cal(pore_length,sample_start_position) for _ in range(500)]
print(f"sampled energy {sampled_energy}")

#initial_energies = [0,1,10,100,700]  # in eV
initial_energies = sampled_energy
Q_m = 9.58e7
ion_mass = 938.272e6  # Mass of proton in eV/c^2
#ions = [Q_m,Q_m/4,Q_m/8,Q_m/16,Q_m/32,Q_m/64]
#ion_mass_in_ev = [ion_mass,ion_mass*4,ion_mass*8,ion_mass*16,ion_mass*32,ion_mass*64]
ions = [Q_m/18]
ion_mass_in_ev = [ion_mass*18]
norm_ion_mass = [mass / 938.272e6 for mass in ion_mass_in_ev]
voltages = [700]  # in V
times_by_energy_voltage = {energy: {voltage: [] for voltage in voltages} for energy in initial_energies}


plt.figure(figsize=(16, 8))
for mass, charge_mass_ratio in zip(ion_mass_in_ev, ions):
    for voltage in voltages:
        for initial_energy in initial_energies:
            t = 0
            xi, yi, ti = electron_trajectory(target_distance, charge_mass_ratio, mass, initial_energy, 8, voltage)
            if initial_energy != intital_energy_ion:
                # where 0 is the angle of the electric field in the pore and 8 is reffering 
                # to the angle of field in the mcp gap relative to the pore

                x, y, t = electron_trajectory(pore_length/(750/initial_energy), charge_mass_ratio, mass, intital_energy_ion,
                                            0, initial_energy)
                
                x_rotated, y_rotated = rotate_coordinates(x, y, 8)
                # Convert to mm
                x_rotated = np.array(x_rotated-x_rotated[-1]) * 1e3
                y_rotated = np.array(y_rotated-pore_length/(750/initial_energy)) * 1e3

                plt.plot(x_rotated, y_rotated, label=f'{mass / 1e6:.1f} MeV/c^2, t_o_f in pore: {t[-1] / 1e-9:.4f},'
                                                    f' int_en: {initial_energy} ev')
                print(f"total time {t[-1]*1e9} ns")
                t = t[-1]
            if initial_energy == intital_energy_ion: 
                xi, yi, ti = electron_trajectory(target_distance, charge_mass_ratio, mass, initial_energy, 0, voltage)
            xi = np.array(xi) * 1e3
            yi = np.array(yi) * 1e3
            total_t = ti[-1] + t
            if t == 0:
                print("this was from rest, top of the pore")
            else:
                print(f"this was from the bottom of the pore, then took {ti[-1]*1e9} ns to get to the top")

            print(f"total time {total_t*1e9} ns")
            plt.plot(xi, yi, label=f'{mass / 1e6:.1f} MeV/c^2, t_o_f: {ti[-1] / 1e-9:.4f}'
                                f', int_en: {initial_energy} ev')
            times_by_energy_voltage[initial_energy][voltage].append(total_t * 1e9)


plt.xlabel('x position (mm)')
plt.ylabel('y position (mm)')


# Title and Legend
plt.title('Electron Trajectories for Different Ion Masses,\n coming out at a 8 degree angle with [0,1,10,100,1000] ev')
#plt.legend()
plt.figure(figsize=(16, 8))
colors = ['b', 'g', 'r','c','m','y','k']
#for energy, color in zip(initial_energies, colors):
   # for idx, voltage in enumerate(voltages):
   #     plt.plot(norm_ion_mass, times_by_energy_voltage[energy][voltage], 'o-', color=colors[idx])
        
#for idx, voltage in enumerate(voltages):
#    min_times = np.min([times_by_energy_voltage[energy][voltage] for energy in initial_energies], axis=0)
#    max_times = np.max([times_by_energy_voltage[energy][voltage] for energy in initial_energies], axis=0)
#    plt.fill_between(norm_ion_mass, min_times, max_times, color=colors[idx], alpha=0.2, label=f'{voltage} V')




plt.xlabel('Mass number')
plt.ylabel('Time [ns]')
plt.title('Time of Flight spread vs. Mass Number at different voltages across the cathode MCP gap ')
#plt.legend()
plt.figure(figsize=(16, 8))
ion_mass_index = 0  # Index of the ion mass you want to plot
times_for_one_ion = []
for energy in initial_energies:
    times = times_by_energy_voltage[energy][voltage][ion_mass_index] 
    times_for_one_ion.append(times)
   
print(f"len {len(times_for_one_ion)}")
plt.hist(times_for_one_ion, bins=50)
times_for_one_ion_array = np.array(times_for_one_ion)  

print(f"type {type(times_for_one_ion_array)}")

plt.xlabel('Time [ns]')
plt.ylabel('Counts')
plt.title(f'Histogram of Time for Ion Mass 1u at  0-1000ev (10ev) steps')
#plt.yscale('log')
print(f"times {times_for_one_ion_array}")
plt.show()

np.savetxt("./data_from_ion_pot/times_by_energy_voltage_ar.csv", times_for_one_ion_array, delimiter=",")





#paper_data = pd.read_csv('Ion_feedback_data/paper/mass_time_table.txt', sep=',', names=['Source', 'Ion',
#                                        'Mass [kg]', 'Calculated time [μs]', 'Corresponding after-pulse group',
#                                        'Average measured time [μs]'], skiprows=1)

#calculated_times = np.array(paper_data['Calculated time [μs]'])

# Perform the division of the DataFrame by the last value in the 'Calculated time [μs]' column
#paper_data_normalized = paper_data.copy()  # Create a copy to avoid modifying the original data
#paper_data_normalized['Calculated time [μs]'] = paper_data['Calculated time [μs]'] / calculated_times[1]

#print(paper_data_normalized['Calculated time [μs]'])

#sampling_rate = 3e9  # 3 GS/s (3 billion samples per second)


# Step 2: Generate some sample data (for example, a sine wave)


# Step 3: Define the bin size based on the time interval
#time_per_sample = 1 / sampling_rate  # Time between each sample (333.33 ps for 3 GS/s)
#bin_size = 10e-9  # Each bin represents 10 nanoseconds

# Convert bin size from time to number of samples per bin
#num_samples_per_bin = int(bin_size * sampling_rate)


# os = pd.read_csv("time_differences.csv", header=None)
# os_1 = pd.read_csv("time_differences_1.csv", header=None)
# plt.figure(figsize=(16, 8))
# plt.hist(os, bins=len(os)//num_samples_per_bin ,alpha=0.5, label="1540V ions")
# plt.hist(os_1, bins=len(os)//num_samples_per_bin, alpha=0.5, label="1780V ionS")
# #plt.xlim(0, 20)
# plt.legend()
# plt.xlabel('Time difference [ns]')
# plt.ylabel('Frequency')
# plt.show()
