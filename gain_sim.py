import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
from scipy.linalg import norm

from energy_dis import accept_reject_v, sey_coefficient, inverse_cdf_output, universal_yield_curve,sey_coefficient_guest
from stepping import electron_yield, solve_for_intercept_time, step_position, step_energy, step_velocity
from consin_dis import cosine_dis, calculate_theta_cylinder

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

def generate_electron_yield(mean_yield: float) -> int:
    # Base number of electrons emitted (integer part)
    base_electrons = int(mean_yield)
 
    # Probability of emitting an additional electron (decimal part)
    additional_prob = mean_yield - base_electrons

    # Roll a random number between 0 and 1
    random_roll = np.random.random()

    # If the roll is less than the probability, emit an extra electron
    if random_roll < additional_prob:
        total_electrons = base_electrons + 1
    else:
        total_electrons = base_electrons

    return total_electrons



def run_simulation(r, v, starting_energy, starting_angle, single_run=True):
    # Initialize variables that were previously global
    first_bounce = True
    max_electrons = 1000000
    @dataclass
    class Electron:
        position: List[float]
        impact_energy : float
        angle: float
        bounce_id: int
        initiating_electron: Optional['Electron'] = None
        secondaries: List['Electron'] = field(default_factory=list)
        end_point: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
        yield_term: float = 9.0
        emitted_energy: float = 0
        # To store the end point

        def parent(self):
            return self.initiating_electron

        def gen_secondaries(self, end_point, impact_energy,angle):
            # Generate secondary electrons
            secondaries = []
        
            poisson_mean = sey_coefficient(impact_energy, angle)
            #poisson_mean = sey_coefficient_guest(impact_energy, angle )
            
            #poisson_mean = int(np.round(poisson_mean))
            #poisson_mean = np.random.poisson(poisson_mean, 1)
            poisson_mean = generate_electron_yield(poisson_mean)

            #poisson_mean = poisson_mean[0]
            electrons_yield = poisson_mean
            self.yield_term = poisson_mean

            total_assigned_energy = 0
            available_energy = self.impact_energy  # Track remaining energy
            delta = poisson_mean  
            
            

            if end_point[2] > end_position_threshold:
                electrons_yield = 0
                self.yield_term = 0
            
            for _ in range(electrons_yield):
                initial_energy_distribution = inverse_cdf_output(1000, available_energy, T, delta)  # if electrons_yield > 0 else 0
                int_energy_from_emitted_electron = random.choice(initial_energy_distribution)

                # Ensure sum of energies does not exceed impact_energy
                if total_assigned_energy + int_energy_from_emitted_electron > self.impact_energy:
                    print(f"sum of energies is higher then impact_energy, electrons {delta}")
                    break  # Stop generating secondaries
                
               
                available_energy -= int_energy_from_emitted_electron
                
                total_assigned_energy += int_energy_from_emitted_electron
                if available_energy < 0:
                    print(f"error {available_energy} ")

                #emitted_energy = self.emitted_energy 

                secondary = Electron(
                    position=end_point,
                    impact_energy = int_energy_from_emitted_electron,
                    #impact_energy=self.impact_energy,  # Example impact_energy reduction
                    angle=self.angle,  # Example angle change
                    bounce_id=self.bounce_id + 1,
                    initiating_electron=self

                    
                )
                secondaries.append(secondary)
            
            return secondaries


    def track_electron(electron, r, v, m=511e3):
        """
        Track an electron and determine its end point and impact energy
        """

        c = scipy.constants.speed_of_light * 1e-6
        x0 = electron.position  # Initial position  
        d = r*2*71.831
        #d = r*2*46
        E = v * (c ** 2) / (d * m) 
        field_orientation = np.array([0., 0., 1])
        #field_orientation = np.array([np.sin(0.13962634), 0, np.cos(0.13962634)])
        a0 = E * field_orientation

        if electron.bounce_id == 0:
            starting_electron_energy = electron.impact_energy
            angle_of_pore = electron.angle 
            print(f"angle of pore {np.rad2deg(angle_of_pore)}")
            start_orientation = np.array([np.sin(angle_of_pore), 0., np.cos(angle_of_pore)])
            v = np.sqrt(2 * starting_electron_energy / m) * c
            
            v0 = v * start_orientation
           
            t = solve_for_intercept_time(x0, v0, a0, r)
            
            hit_position = step_position(x0, v0, a0, t)
            
            impact_energy = step_energy(v0, a0, t, m)
            print(f"impact energy  {impact_energy}")
            
            end_velocity = step_velocity(v0, a0, t)
            
            impact_angle, impact_angle_d = calculate_theta_cylinder(end_velocity, hit_position, r)
            
        
        else:
            int_energy_from_emmited_electron = electron.impact_energy
            #delta = 1
            #initial_energy_distribution = inverse_cdf_output(10, E_im, T, delta)
            #int_energy_from_emmited_electron = random.choice(initial_energy_distribution)
            
            v_initial = np.sqrt((int_energy_from_emmited_electron / 511e3) * 2 * (c ** 2))
            vel_a, theta = cosine_dis(electron.position , r)
            electron.angle = theta
            v1 = np.array([vel_a[0][0] * v_initial, vel_a[0][1] * v_initial, vel_a[0][2] * v_initial])
            ti = solve_for_intercept_time(electron.position , v1, a0, r)
            
            hit_position = step_position(electron.position, v1, a0, ti)
            
            impact_energy = step_energy(v1, a0, ti, m)
            end_velocity = step_velocity(v1, a0, ti)
            impact_angle, impact_angle_d = calculate_theta_cylinder(end_velocity, hit_position, r)
           


        # Track electron and determine its end point
        # Simplified example: calculate the end point based on start position and angle
        end_point = hit_position
        electron.end_point = end_point
        electron.impact_energy = impact_energy
        electron.angle = impact_angle
        return  end_point,impact_energy, impact_angle


    def generate_secondary(end_point, impact_energy, angle, electron):
        """
        Generate secondary electrons based on the end point and impact energy of the parent electron
        """
        return electron.gen_secondaries(end_point, impact_energy, angle)


    def track_electron_and_generate_secondaries(electron,r,v):
        """
        Track an electron and generate secondary electron
        """
        end_point,impact_energy, angle = track_electron(electron,r,v)
        
        secondaries = generate_secondary(end_point, impact_energy, angle, electron)
        electron.secondaries = secondaries
        return secondaries


    def gen_electron():
        # Generate the first electron
        return Electron(position=[0, 0, 0], impact_energy=starting_energy, angle=np.deg2rad(starting_angle), bounce_id=0, yield_term=1.0)

    electrons = []
    if first_bounce:
        first_electron = gen_electron()
        electrons.append(first_electron)

    all_electrons = []

    total_electron_count = 0
    end_position_count = 0
    r= r
    end_position_threshold = r*2*71.83 #38.3
    T = 2.5  # eV, assumed temperature
    
    start_positions = []
    end_positions = []
    energies = []
    energies_out = []
    out_of_pore = []
    in_pore = []
    electron_yield_values = []


    while electrons:
        # takes the last electron in the list from electrons
        e = electrons.pop()
        bounce_id = e.bounce_id
    
        new_electrons = track_electron_and_generate_secondaries(e,r,v)
        total_electron_count += len(new_electrons)

        
        if total_electron_count >= max_electrons:
            break  # Stop the simulation if the maximum number of electrons is exceeded

        if e.end_point[2] > end_position_threshold:
            out_of_pore.append(e)
            energies_out.append(e.impact_energy)
        if e.end_point[2] < end_position_threshold:
            in_pore.append(e)
            energies.append(e.impact_energy)
            electron_yield_values.append(e.yield_term)

        all_electrons.extend(new_electrons)
        electrons.extend(new_electrons)
        start_positions.append(e.position)
        end_positions.append(e.end_point)
        
        
    print(f"Total electrons generated: {total_electron_count}")
    #print(f"yields {np.mean(electron_yield_values)}")

    # Plotting impact_energy vs. bounces
    electron_yield_values_list = []
    for i in energies:
        poisson_mean = sey_coefficient(i, 0)
        #poisson_mean = sey_coefficient_guest(i,0)
        #poisson_mean = int(np.round(poisson_mean))
        
        poisson_mean = np.random.poisson(poisson_mean, 1)
        poisson_mean = poisson_mean[0]
        electron_yield_values_list.append(poisson_mean)
    #print(f"electron_yield_values_list {electron_yield_values_list}")
    

    #plt.hist(electron_yield_values_list, bins=100, alpha=0.75, range=(0,100), label='Electron Energy')

    bounce_counts = {}
    for electron in all_electrons:
        bounce_id = electron.bounce_id
        
        if bounce_id in bounce_counts:
            bounce_counts[bounce_id] += 1
        else:
            bounce_counts[bounce_id] = 1

    start_positions = np.array(start_positions)
    end_positions = np.array(end_positions)

    end_positions_np = np.array(end_positions)
    count_above_threshold = np.sum(end_positions_np[:, 2] > end_position_threshold)

    print(f"Number of end positions above the threshold {end_position_threshold}: {count_above_threshold}")

    print((bounce_counts))
    if single_run == True:

        plt.figure(figsize=(12, 6))
        plt.scatter(energies,electron_yield_values_list, color='blue')
        plt.xlabel('Electron Energy')
        plt.ylabel('Electron Yield')
        plt.title('Electron Yield vs. Electron Energy')
        plt.figure(figsize=(12, 6))
        plt.hist(electron_yield_values_list, bins=100, alpha=0.75, range=(0,4), label='Electron Energy')
        plt.xlabel('Electron Yield')
        plt.ylabel('Number of Electrons')
        plt.title(f'Electron Yield Distribution mean = {np.mean(electron_yield_values_list)}')


        # Prepare data for plotting
        bounces = list(bounce_counts.keys())
        electron_counts = list(bounce_counts.values())
        fig = plt.figure(figsize=(16, 8))  # Increase the size of the figure
        ax = fig.add_subplot(111, projection='3d')

        # Plot start positions
        ax.scatter(start_positions[:, 0], start_positions[:, 2], start_positions[:, 1], color='blue', label='Start Positions')
        # Plot end positions
        ax.scatter(end_positions[:, 0], end_positions[:, 2], end_positions[:, 1], color='red', label='End Positions')
        ax.scatter(end_positions[0, 0], end_positions[0, 2], end_positions[0, 1], color='green', label='fist bounce')

        ax.set_xlabel('X')
        ax.set_ylabel('z')
        ax.set_zlabel('y')
        ax.set_ylim(0, end_position_threshold)
        ax.set_title('Electron Start and End Positions')
        ax.legend()

        # Cylinder parameters
        origin = np.array([0, 0, 0])
        p0 = np.array([0, 1e-3, 0])
        p1 = np.array([0, end_position_threshold, 0])
        R = r  # Radius of the cylinder
        print(f"end position threshold {end_position_threshold}")

        # Vector in the direction of the axis
        v = p1 - p0
        mag = norm(v)
        v = v / mag

        # Two perpendicular vectors to the axis vector
        not_v = np.array([1, 0, 0])
        if np.allclose(v, not_v):
            not_v = np.array([0, 1, 0])

        n1 = np.cross(v, not_v)
        n1 /= norm(n1)
        n2 = np.cross(v, n1)

        # Surface ranges
        theta = np.linspace(0, 2 * np.pi, 100)
        z = np.linspace(0, mag, 100)
        theta, z = np.meshgrid(theta, z)

        # Generate coordinates for surface
        X = p0[0] + R * np.sin(theta) * n1[0] + R * np.cos(theta) * n2[0] + v[0] * z
        Y = p0[1] + R * np.sin(theta) * n1[1] + R * np.cos(theta) * n2[1] + v[1] * z
        Z = p0[2] + R * np.sin(theta) * n1[2] + R * np.cos(theta) * n2[2] + v[2] * z

        colors = np.where(((Y >= 0) & (Y <= 1e-2)) | ((Y >= end_position_threshold-0.01) & (Y <= end_position_threshold)), 'red', 'aqua')

        # Plot the surface
        ax.plot_surface(X, Y, Z, facecolors=colors, alpha=0.02)
        ax.set_box_aspect([1, 5, 1])  # Aspect ratio is 1:1:2

        print(f"electron counts {len(end_positions)}")
        plt.figure(figsize=(12, 6))

        plt.bar(bounces, electron_counts, edgecolor='w')
        plt.xlabel('Bounce Count')
        plt.ylabel('Number of Electrons')
        plt.title('Total Number of Electrons per Bounce')
        plt.grid(True)

        plt.figure(figsize=(12, 6))
        #plt.scatter(end_positions[:, 0], end_positions[:, 2], color='blue', label='Start Positions')
        print(f"end positions {end_positions[:,2]}")
        bins = np.arange(0, end_position_threshold + 0.005, 0.005)
        
        hist_vals, bin_edges, _  = plt.hist(end_positions[:,2], bins=bins, color='red', label='End Positions')
        plt.hist(end_positions[:,2], bins=bins, color='red', label='End Positions')
        plt.scatter(bin_edges[:-1]+0.005/2, hist_vals, color='green', label='End Positions')    
        plt.xlim(0,end_position_threshold)
        plt.xlabel('Z Position mm')
        plt.ylabel('Number of Electrons')
        plt.title('Electron hit Position Distribution in the Z direction')

        data_to_save = np.column_stack((bin_edges[:-1]+0.005/2, hist_vals))
        #np.savetxt("end_positions_2.txt", data_to_save,  delimiter=",")
        #existing_data = np.loadtxt("C:/Users/lexda/PycharmProjects/electron_pore_trajectories/end_positions_2.txt", delimiter=",")
        # combined_data = np.vstack((existing_data, data_to_save))
        # combined_data = combined_data[np.argsort(combined_data[:, 0])]
        # unique_bin_edges, indices = np.unique(combined_data[:, 0], return_index=True)
        # summed_values = np.add.reduceat(combined_data[:, 1], indices)
        # updated_data = np.column_stack((unique_bin_edges, summed_values))
        #np.savetxt("end_positions_2.txt", updated_data, delimiter=",", comments='')


        plt.figure(figsize=(12, 6))
        plt.hist(energies, bins=100, alpha=0.75, range=(0,1000), label='Electron Energy')
        plt.hist(energies_out, bins=100, alpha=0.75,range=(0,1000), label='Electron Energy out of pore')
        print(f"energies {np.max(energies)}")
        try:
            print(f"energies_out {np.max(energies_out)}")
        except Exception as e:
            print("no electrons out of pore")
        plt.range = (0, 2000)
        plt.xlabel('Electron Energy')
        plt.ylabel('Number of Electrons')
        plt.title('Electron Energy Distribution')
        plt.legend()
    
    return  count_above_threshold

gain_spread_mean_list = []

#range_of_r = np.linspace(0.0005, 0.002, num=10)
#range_of_r =  [0.001665,0.00333,0.004,0.001,0.005]
range_of_r = [0.0085]

#range_of_v = [1200,1000,800,600]
range_of_v = [1000]
angles = [10] #list(range(2, 40, 5))
number_of_runs = 1

for v in range_of_v:
    for r in range_of_r:
        
        starting_energy =  200
        starting_angle = angles[0]
        gain_spread = []
        for j in range(number_of_runs):
            count_above_threshold= run_simulation(r,v,starting_energy,starting_angle,single_run=True)
            gain_spread.append(count_above_threshold)
            
        plt.hist(gain_spread, bins=50, alpha=0.75, label='Electron Energy',log=True)
        plt.xlabel('Gain (Number of Electrons leaving one pore)')
        plt.ylabel('counts')
        gain_spread_mean = np.mean(gain_spread)
        gain_spread_mean_list.append((v, r, gain_spread_mean))
    
print(f"gain spread mean {gain_spread_mean}")

print(f"gain spread {gain_spread_mean_list}")
single_run = False
if single_run == False:  
    plt.figure(figsize=(12, 6))

    for v in range_of_v:
        r_values = np.array( [r for v_val, r, gain in gain_spread_mean_list if v_val == v])
        gain_values = [gain for v_val, r, gain in gain_spread_mean_list if v_val == v]
        plt.scatter((0.00165*2*71.83)/(r_values*2), gain_values, label=f'v = {v:.2f}')

    plt.ylabel('Number of Electrons leaving the pore')
    plt.xlabel('L/D')
    plt.legend()
    plt.yscale('log')
plt.show()



