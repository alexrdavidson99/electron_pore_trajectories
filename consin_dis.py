"""
sampling from a cosine distribution
https://www.particleincell.com/2015/cosine-distribution/

@author: lubos brieda
"""
from random import random
import math
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def cosine_dis(xi,r):
    """
    Sample random points from a cosine distribution.
    """
    angle_dis = []
    vel_a = []
    
    # surface properties
    tang1 = [-xi[1]/r, xi[0]/r, 0]
    tang2 = [0, 0, -1]
    norm = [-xi[0]/r, -xi[1]/r, 0]

    # sample random points
    for it in range(0, 1):
        sin_theta = math.sqrt(random())
        cos_theta = math.sqrt(1 - sin_theta * sin_theta)

        # random in plane angle
        psi = random() * 2 * math.pi;

        # three vector components
        a = sin_theta * math.cos(psi)
        b = sin_theta * math.sin(psi)
        c = cos_theta

        angle_dis.append(a)

        # multiply by corresponding directions
        v1 = [a * n for n in tang1]
        v2 = [b * n for n in tang2]
        v3 = [c * n for n in norm]

        # add up to get velocity, vel=v1+v2+v3

        vel = []
        for i in range(0, 3):
            vel.append(v1[i] + v2[i] + v3[i])

        vel_a.append(vel)

        #print(np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2))
        theta = math.acos(cos_theta) * 180 / math.pi;



    return vel_a, theta



import numpy as np

def calculate_theta_cylinder(velocity, xi, r):
    
    x = xi[0]
    y = xi[1]
    z = xi[2]
    # Define the normal vector to the x-z plane (along the y-axis)
    normal_vector = np.array( [x / r, y / r, 0])
    
    # Normalize the velocity vector
    velocity_magnitude = np.linalg.norm(velocity)
     
    velocity_normalized = velocity / velocity_magnitude
    
    # Normalize the normal vector (although it's already a unit vector, just for completeness)
    normal_magnitude = np.linalg.norm(normal_vector)
    
    normal_normalized = normal_vector / normal_magnitude
    
    # Compute the dot product between the velocity vector and the normal vector
    dot_product = np.dot(velocity_normalized, normal_normalized)
    
    # Ensure the dot product stays within the valid range [-1, 1] for arccos
    dot_product = np.clip(dot_product, -1.0, 1.0)
    
    # Calculate the angle in radians
    theta = np.arccos(dot_product)
    
    # Convert the angle to degrees
    theta_degrees = np.degrees(theta)
    
    return theta, theta_degrees
   



plot = False

if plot:

    xi = [0.2724217, 0., 0.23779032]
    r = 0.02724217
    v_bin = []
    
    for i in range(0, 10000):
        vel_a, theta = cosine_dis(xi, r)
        

        v_en = 1
        v1 = np.array([vel_a[0][0]*v_en, vel_a[0][1]*v_en, vel_a[0][2]*v_en])
        v_bin.append(v1)
    
    v_bin = np.array(v_bin)
    fig = plt.figure()
    ax2 = fig.add_subplot(212, projection='3d')
    ax2.scatter(v_bin[:, 0],v_bin[:, 1],v_bin[:, 2],c='r',s=0.005);
    ax2.set_xlabel('angle')
    ax2.set_ylabel('normalized count')
    plt.show()


