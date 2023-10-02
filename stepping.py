import numpy
import numpy as np
import scipy
import matplotlib.pyplot as plt


def solve_for_intercept_time( x0, v0, acc, radius ):
    '''
    Solve for intercept on the cylinder
    ''' 
    _coeff = [ 0.25*(acc[0]**2 + acc[1]**2),  # t^4
               ( acc[0]*v0[0] + acc[1]*v0[1]),  # t^3
               ( v0[0]**2 + v0[1]**2 + x0[0]*acc[0] + x0[1]*acc[1] ),  # t^2
               2.0*( x0[0]*v0[0] + x0[1]*v0[1] ),  # t
               ( x0[0]**2 + x0[1]**2 - radius**2 ) ]  # constant

    _roots = numpy.roots( _coeff )

    _roots = _roots[ numpy.where( numpy.imag( _roots ) == 0, True, False ) ]
    _roots = _roots[(_roots > 0) & (_roots > 1e-15)]


    return numpy.real( numpy.min( _roots ) )


def step_position( x0, v0, acc, time ):
    '''
    Position after time step
    '''
    return x0 + v0*time + 0.5*acc*time**2

def step_energy( v0, acc, time ):
    '''
    Energy after time step starting from emission velocity
    '''
    m = 511e3 # in eV
    c = scipy.constants.speed_of_light*1e-6 # in mm/ns
    return 0.5*m*(numpy.dot(v0 + acc*time, v0 + acc*time))/(c**2)


def solve_for_exit_time(x0,v0, acc, depth):
    '''
    Exit time, starting from the emission position
    '''
    return (-v0[2] + numpy.sqrt(v0[2]**2 - 2*acc[2]*(x0[2] - depth)))/(acc[2])


if __name__ == '__main__':
    
    c = scipy.constants.speed_of_light*1e-6 # in mm/ns
    V = 1000.   # in V
    d = 2e-1    # pore depth in mm
    r = 2e-3    # pore radius in mm
    m = 511e3   # in eV
    E = V*(c**2)/(d*m) # electric field acceration in mm/ns^2
    e = 10      # in eV
    v = numpy.sqrt( 2*e/m )*c # velocity in mm/ns

    #
    orientation = numpy.array([ numpy.sin( 0.2 ), 0., numpy.cos( 0.2 ) ] )
    
    x0 = numpy.array( [ 0, 0, 0 ] )
    v0 = v*orientation

    a0 = E*orientation


    time = []
    position_y = []
    position_x = []
    energy = []
    velocity = []
    total = 0


    # Step particle in field until it reaches boundary or exits
    # inital step
    t = solve_for_intercept_time( x0, v0, a0, r )
    print(t)
    x1 = step_position(x0, v0, a0, t)
    print(x1, x1[0] ** 2 + x1[2] ** 2)
    energy1 = step_energy(v0, a0, t)

    #inital secondaries
    v1 = np.array([-np.sqrt((3.95/(511e3))*2*(c**2)), 0, 0])

    for i in range(3):
        ti = solve_for_intercept_time(x1, v1, a0, r)
        xi = step_position(x1, v1, a0, ti)
        for j in range(100):

            xj = step_position(x1, v1, a0, (ti/100)*(j))
            position_x.append(xj[0])
            position_y.append(xj[2])



        print(xi, xi[1] ** 2 + xi[2] ** 2)
        energy2 = step_energy(v1, a0, t)
        print(energy2)
        x1 = xi
        #v1 = -v1
        total += ti

        #position.append(xi[2])
    print(position_x)
    print(position_y)
    print(total)
    plt.plot(position_x,position_y)
    #plt.scatter(range(len(time)),time)
    #plt.scatter(range(len(position_y)),position_y)
    plt.title(f"Electron position over time, with initial energy of 3.95eV where r is {r} and d is {d} mm")
    plt.xlabel("x position (mm)")
    plt.ylabel("y position (mm)")
    plt.hlines(d, -r, r, colors='k', linestyles='solid')
    plt.show()



    #print(total)
    #print(position)
    #print(position[5]-position[4])
    #print(position[7] - position[6])
















