import numpy
import numpy as np
import scipy
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import photutils
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from scipy.stats import cosine
from scipy.linalg import norm
from consin_dis import cosine_dis
import functions_for_saving as ffs



class Electron:
    def __init__(self, initial_position, initial_velocity):
        self.position_x = [initial_position[0]]
        self.position_y = [initial_position[2]]
        self.position_z = [initial_position[1]]
        self.velocity = initial_velocity

    def generate_new_electrons(self, e_energy):
        # Implement your logic to generate new electrons based on the given energy
        # Example: Generating new electrons based on a yield function
        print("test_yield_energy")

        electrons_yield = electron_yield(e_energy)
        new_electrons = [Electron([self.position_x[-1], self.position_z[-1], self.position_y[-1]], self.velocity) for _
                         in range(electrons_yield)]
        return new_electrons



def p(x):

    # return x/(5**2)*np.exp(-(x**2)/(2*(5**2)))
    """
    This is the pdf for the energy distribution, it is a log normal distribution.
    :param x:
    :return: pdf
    """
    # return (1 / (x * 1.0828 * (2 * np.pi) ** (1 / 2))) * np.exp(-(((np.log(x) - 1.1) ** 2) / (2 * (1.0828 ** 2))))
    return (x / (1 ** 2)) * np.exp(-x / 1)


def accept_reject(N):
    xmin = 0
    xmax = 250
    # pmax = 0.12130
    # pmax = (1 / 1) * np.exp(-1)
    pmax = 0.12544838348768622

    n_accept = 0
    x_list = []
    while n_accept < N:
        t = (xmax - xmin) * np.random.rand() + xmin
        y = np.random.rand()
        if y < p(t) / pmax:
            n_accept += 1
            x_list.append(t)
    return x_list


def solve_for_intercept_time( x0, v0, acc, radius ):
    '''
    Solve for intercept on the cylinder
    '''

    print(x0[0]*v0[0], x0[1]*v0[1])
    print(x0[1])
    _coeff = [ 0.25*(acc[0]**2 + acc[1]**2),  # t^4
               ( acc[0]*v0[0] + acc[1]*v0[1]),  # t^3
               ( v0[0]**2 + v0[1]**2 + x0[0]*acc[0] + x0[1]*acc[1] ),  # t^2
               2.0*( x0[0]*v0[0] + x0[1]*v0[1] ),  # t
               ( x0[0]**2 + x0[1]**2 - radius**2 ) ]  # constant
    print(_coeff)
    _roots = numpy.roots( _coeff )
    print(_roots)

    _roots = _roots[ numpy.where( numpy.imag( _roots ) == 0, True, False ) ]

    _roots = _roots[(_roots > 0) & (_roots > 1e-16)]

    print(numpy.real( numpy.min( _roots ) ))

    return numpy.real( numpy.min( _roots ) )


def step_position( x0, v0, acc, time ):
    '''
    Position after time step
    '''
    return x0 + v0*time + 0.5*acc*time**2

def step_energy( v0, acc, time, m ):
    '''
    Energy after time step starting from emission velocity
    '''
    c = scipy.constants.speed_of_light*1e-6 # in mm/ns
    return 0.5*m*(numpy.dot(v0 + acc*time, v0 + acc*time))/(c**2)

def electron_yield(energy2):
    # poisson_mean = 4.25 * (energy1 / 600) * np.exp(1 - (energy1 / 600))
    #poisson_mean = 4.25 * ((energy1 / 600) * np.exp(1 - energy1 / 600)) ** 0.56
    poisson_mean = 4.25 * (energy2 / 600) * np.exp(1 - (energy2 / 600))
    print(f"poisson mean {poisson_mean}")
    pos_0 = np.random.poisson(poisson_mean, 100)
    electrons_yield = np.random.choice(pos_0)
    print(f"yield {electrons_yield}")
    if electrons_yield < 1:
        print(f"no more electrons, yield  {electrons_yield}")
    s_yield.append(electrons_yield)
    return  electrons_yield


def solve_for_exit_time(x0,v0, acc, depth):
    '''
    Exit time, starting from the emission position
    '''

    return (-v0[2] + numpy.sqrt(v0[2]**2 - 2*acc[2]*(x0[2] - depth)))/(acc[2])


if __name__ == '__main__':
    plt.ioff()
    all_n_of_collisions = []
    x_end = []
    y_end = []
    energy_overall = []
    total_yield =[]
    accept_reject_energy = accept_reject(100000)
    one_run = False # set to true if you want to run one simulation
    e_field_in_z = True # set to true if you want the e-field to be in the z direction
    message_printed = False # used to print the message only once

    for _ in range(3):

        c = scipy.constants.speed_of_light*1e-6 # in mm/ns
        V = 1000.   # electrods potential in V
        d = 2e-1    # pore depth in mm
        l_over_d = 40
        r = (d/l_over_d)*0.5 
       
      
        #r = 1.666666e-3    # pore radius in mm
        m = 511e3   # in eV
        #m = 195303.27e6
        E = V*(c**2)/(d*m)  # electric field acceleration in mm/ns^2
        e = 200  # in eV
        v = numpy.sqrt(2*e/m)*c  # velocity in mm/ns
        print(1.6e-19 * 2 / 2 * r)
        print(f"l/D = {(d/(2*r)):.2f} length to diameter ratio")
        orientation = numpy.array([numpy.sin(0.13962634), 0., numpy.cos(0.13962634)])

        #x0 = numpy.array([-r, 0, 0])
        x0 = numpy.array([0, 0, 0])
        v0 = v*orientation

        if e_field_in_z is True:
            orientation = numpy.array([0., 0., 1])
        a0 = E*orientation

        time = []
        n_of_collisions = []
        x_hit = []
        y_hit = []
        z_pos = []
        energy = []
        s_yield = []
        error = []
        angle = []
        velocity = []
        total = 0

        # Step particle in field until it reaches boundary or exits
        # initial step
        t = solve_for_intercept_time(x0, v0, a0, r)
        x1 = step_position(x0, v0, a0, t)
        print(x1)

        energy1 = step_energy(v0, a0, t, m)
        energy.append(energy1)

        electrons_yield = electron_yield(energy1)

        x_hit.append(x1[0])
        y_hit.append(x1[1])
        z_pos.append(x1[2])
        time.append(t)

        # initial secondaries
        e_p = random.choice(accept_reject_energy)
        print(f"energy from accept_reject {e_p}")

        vel_a,theta = cosine_dis(x1, r)
        angle.append(theta)
        angle_radians = np.arccos(x1[0] / r)
        print(f" angle need for rotation {angle_radians * (180 / np.pi)}")
        v_en = np.sqrt((e_p / 511e3) * 2 * (c ** 2))
        v1 = np.array([vel_a[0][0] * v_en, vel_a[0][1] * v_en, vel_a[0][2] * v_en])
        print(f" v1 before rotation {v1[0]} 2 {v1[1]} 3 {v1[2]}")
        num_electrons = electrons_yield

        electrons = [Electron(x0, v0) for _ in range(num_electrons)]
        initial_electrons = electrons.copy()

        # number of collisions with the wall  - 1
        for current_electron in initial_electrons:
            x1 = x0
            v1 = v0

            current_electron.position_x = [x0[0]]
            current_electron.position_y = [x0[2]]
            current_electron.position_z = [x0[1]]
            current_electron.velocity = v0

            for i in range(206):     #206
                ti = solve_for_intercept_time(x1, v1, a0, r)
                print(f" time to hit {ti}")
                xi = step_position(x1, v1, a0, ti)

                for j in range(100):
                    xj = step_position(x1, v1, a0, (ti/100)*j)
                    current_electron.position_x.append(xj[0])
                    current_electron.position_y.append(xj[2])
                    current_electron.position_z.append(xj[1])

                    #out of pore condition
                    if xj[2] > d and not message_printed:
                        print(f"out of pore in {i - 1} collisions")
                        message_printed = True
                        #n_of_collisions.append(i + 1)
                        x_end.append(xj[0])
                        y_end.append(xj[1])

                print(xi, xi[0] ** 2 + xi[1] ** 2)
                energy2 = step_energy(v1, a0, ti, m)
                energy.append(energy2)

                electrons_yield = electron_yield(energy2)

                for _ in range(electrons_yield):
                    e_p = random.choice(accept_reject_energy)
                    vel_a, theta = cosine_dis(xi, r)
                    v_en = np.sqrt((e_p / 511e3) * 2 * (c ** 2))
                    new_velocity = np.array([vel_a[0][0] * v_en, vel_a[0][1] * v_en, vel_a[0][2] * v_en])

                    print(f"This is the new {electrons_yield} {new_velocity} {xi}")
                    new_electron = Electron(xi, new_velocity)
                    electrons.append(new_electron)

                count = 0
                while count < 1000000:
                    e_p = random.choice(accept_reject_energy)
                    if energy2 > e_p:
                        print(f"right energy {energy2} {e_p}")
                        break
                    if count > 999999:
                        e_p = 100
                        print("never got a lower energy")
                        break
                    count += 1

                x_hit.append(xi[0])
                y_hit.append(xi[1])
                z_pos.append(xi[2])

                print(f" end pos {xi} mm")
                print(f" v1 before rotation {v1[0]} 2 {v1[1]} 3 {v1[2]}")

                vel_a, theta = cosine_dis(xi,r)
                angle_radians = np.arccos(xi[0] / r)
                print(f" angle need for rotation {angle_radians*(180/np.pi)}")

                v_en = np.sqrt((e_p/ 511e3) * 2 * (c ** 2))
                v1 = np.array([vel_a[0][0]*v_en, vel_a[0][1]*v_en, vel_a[0][2]*v_en])

                angle.append(np.arccos((v1[0]/np.sqrt(v1[0]**2+v1[1]**2))*(-xi[0])/r + (v1[1]/np.sqrt(v1[0]**2+v1[1]**2))*(-xi[1])/r))

                print(f" v1 after rotation {v1}")
                error.append(v1[0]*-xi[0]/r + v1[1]*-xi[1]/r)
                print(f" is bigger then 0 {v1[0]*-xi[0]/r + v1[1]*-xi[1]/r }")

                x1 = xi
                total += ti
                time.append(total)

                if xi[2] > d:
                    print(f"out of pore in loop {i - 1}")
                    n_of_collisions.append(i)
                    
                    break  # Exit the outer loop when the condition is met
                #if (electron_yield(energy2)) == 0:
                #    print("no more electrons")
                #    break

            message_printed = False # reset the message printed flag for the next electron
        if one_run is True:
            print(n_of_collisions)
            print(total)
            fig = plt.figure(figsize=(10, 12))
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect([1, 2, 1])
            colors = ['red', 'green', 'blue', 'yellow', 'purple', 'orange']  # Add more colors as needed

            for i, new_electron in enumerate(electrons):
                color = colors[i % len(colors)]
                ax.plot(new_electron.position_x, new_electron.position_y, new_electron.position_z, label='Projectile Motion', color=color)
            ax.scatter(x_hit, z_pos, y_hit, color='red', marker='o', label='End Point')

            plt.title(f"Electron position over time, Where E-field is {(V/(d*10**-3)):.2e} r is {r} and d is {d} mm")
            plt.xlabel("x position (mm)")
            plt.ylabel("y position (mm)")

            # # making the cylinder!!!

            # Cylinder parameters
            origin = np.array([0, 0, 0])
            p0 = np.array([0, 1e-3, 0])
            p1 = np.array([0, d, 0])
            R = r  # Radius of the cylinder

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

            colors = np.where(((Y >= 0) & (Y <= 1e-2)) | ((Y >= 1.8e-1) & (Y <= 2e-1)), 'red', 'aqua')

            # Plot the surface
            ax.plot_surface(X, Y, Z, facecolors=colors, alpha=0.02)
            

            plt.figure()
            print(len(time))
            print(len(energy))
            df = pd.DataFrame({"x_hit": x_hit, "y_hit": y_hit, "z_pos": z_pos,"Time": time,  "energy": energy,"angle": angle})
            

            plt.scatter(df.x_hit, df.y_hit, c=df.Time, cmap='viridis', marker='o')
            plt.colorbar(label='Time')



            plt.figure()
            plt.hist(df.angle[1:25], bins=10, alpha=0.7, rwidth=0.85, density=True, color='purple', edgecolor='black')
            pi_ticks = [0, 0.5 * np.pi]
            pi_labels = ['0', 'π/2']
            plt.xticks(pi_ticks, pi_labels)
            plt.title("Cosine dis Histogram of theta ")
            plt.xlabel("Angle (radians)")
            plt.ylabel("probability density")

            plt.figure()
            plt.scatter(df.x_hit, df.z_pos, c=df.energy, cmap='plasma' )#norm=LogNorm(), marker='o')

            plt.hist2d(df.x_hit, df.y_hit, bins=(20, 20), cmap='viridis')

            print(df)
            plt.colorbar(label='Energy (log scale)')
            plt.axhline(y=d, color='red', linestyle='--', label='z_pos = r')

            plt.title("Electron Position with Energy Heat Map")
            plt.xlabel("x_hit (mm)")
            plt.ylabel("z_pos (mm)")

        all_n_of_collisions.append(n_of_collisions)
        energy_overall.append(energy)
        total_yield.append(s_yield)


    flattened_n_of_collisions = [item for sublist in all_n_of_collisions for item in sublist]
    print(f"n... {flattened_n_of_collisions}")
    plt.figure()
    bins_for_collisions = [i + 0.5 for i in range(min(flattened_n_of_collisions) - 1, max(flattened_n_of_collisions) + 1)]
    plt.hist(flattened_n_of_collisions, bins=bins_for_collisions, alpha=0.7, edgecolor='black', align='mid', rwidth=1.0) #, rwidth=0.85)
    flattened_n_of_collisions = np.array(flattened_n_of_collisions)
    
    ffs.append_data_with_header('outputs/Hits_data.txt', f"l/D = {(d / (2 * r)):.2f}", flattened_n_of_collisions)
    plt.title(f"Histogram of n_of_collisions parameters: l/D = {(d/(2*r)):.2f}, E-field = {(V/(d*10**-3)):.2e}, radius = {r}, length = {d}")
    plt.xlabel("Number of Collisions")
    plt.ylabel("Frequency")

    # new plot 
    # Create points on the circle to show phase space
    plt.figure()

    circle_center = (0, 0)
    theta = np.linspace(0, 2 * np.pi, 100)
    circle_x = circle_center[0] + r * np.cos(theta)
    circle_y = circle_center[1] + r * np.sin(theta)

    # Plot the dashed circular line
    plt.plot(circle_x, circle_y, linestyle='dashed', label='pore boundary')
    plt.legend()
    plt.scatter(x_end, y_end, marker='o')
    plt.xlim(-r-r/5, r+r/5)
    plt.ylim(-r-r/5, r+r/5)

    plt.title("pore exit points")
    plt.xlabel("x position (mm)")
    plt.ylabel("y position (mm)")

    # new plot 
    # shows the energy distribution
    plt.figure()
    
    flattened_energy_overall = [item for sublist in energy_overall for item in sublist]
    flattened_yield_overall = [item for sublist in total_yield for item in sublist]
    
    plt.hist(flattened_energy_overall, alpha=0.7, rwidth=0.85,
             edgecolor='black', range=(0, 500), bins=50 , density=True)
    flattened_energy_overall= np.array(flattened_energy_overall)
    ffs.append_data_with_header('outputs/energy_data.txt', f"l/D = {(d / (2 * r)):.2f}", flattened_energy_overall)

    # Calculate percentiles
    percentile_25 = np.percentile(flattened_energy_overall, 25)
    median = np.percentile(flattened_energy_overall, 50)
    percentile_75 = np.percentile(flattened_energy_overall, 75)

    print(f'25th Percentile: {percentile_25}')
    print(f'Median (50th Percentile): {median}')
    print(f'75th Percentile: {percentile_75}')
    print(f"l/D = {(d/(2*r)):.2f} length to diameter ratio")
    print(f"r is {r}")

    plt.title("energy distribution")
    plt.xlabel("energy (eV)")
    plt.ylabel("count")

    # new plot 
    # shows the yeild distribution
    plt.figure()

    
    hist, edges, _ = plt.hist(flattened_yield_overall, alpha=0.7, rwidth=0.85,
             edgecolor='black', density=True, bins=10)

    # Save histogram data to a text file
    hist_data = np.column_stack((edges[:-1], edges[1:], hist))
    np.savetxt('histogram_data_first_hit_CST.txt', hist_data, header='bin_start, bin_end, frequency', fmt='%.4f', delimiter=',')

    plt.title("Histogram of Yield")
    plt.xlabel("Yield")
    plt.ylabel("Frequency")

    print(sum(flattened_yield_overall))

    plt.show()



