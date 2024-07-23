"""
sampling from a cosine distribution
https://www.particleincell.com/2015/cosine-distribution/

@author: lubos brieda
"""
from random import random
import math
import numpy as np



def cosine_dis(xi,r):

    count = []
    angle_dis = []
    vel_a = []
    nbins = 45
    delta_bin = 90 / (nbins)



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




