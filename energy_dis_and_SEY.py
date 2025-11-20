import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import pandas as pd
import random
from scipy.special import gammainc
from consin_dis import calculate_theta_cylinder
from scipy.special import gammainc
from scipy.interpolate import interp1d


def test(impact_energy, angle):
    poisson_mean = sey_coefficient(impact_energy, angle)
    pos_0 = np.random.poisson(poisson_mean, 500)
    return pos_0

def universal_yield_curve(E, Emax=339.79020249876953 , delta_max=3.9283773281400824,sigma= 1.319094301238391):

    delta = delta_max*np.exp((-(np.log(E/Emax))**2)/(2*(sigma**2)))
    return delta

# Defining the PDF with the correct incomplete gamma function parameters
def vaughan_pdf_corrected(E, E0, T, delta):
    # δ(E0, θ0) 
    # Correct incomplete gamma function with (2, E0/T)
    inc_gamma = gammainc(2, E0 / T)
    pdf = delta * (E / T**2) * np.exp(-E / T) / inc_gamma
    return pdf

def p(x,sigmaFit,muFit):
    # return x/(5**2)*np.exp(-(x**2)/(2*(5**2)))
    #sigmaFit = 1.0828
    #muFit = 1.6636
    #return 2.92779 * (x / (7.5 ** 2)) * np.exp(-x / 7.5)
    return (1 / (x * sigmaFit * (2 * np.pi) ** (1 / 2))) * np.exp(-(((np.log(x) - muFit) ** 2) / (2 * (sigmaFit ** 2))))


def accept_reject(N,sigmaFit,muFit):
    xmin = 0
    xmax = 150
    #pmax = 0.12130
    pmax = 0.12544838348768622

    n_accept = 0
    x_list = []
    while n_accept < N:
        t = (xmax - xmin) * np.random.rand() + xmin
        y = np.random.rand()
        if y < p(t,sigmaFit,muFit) / pmax:
            n_accept += 1
            x_list.append(t)
    return x_list


def accept_reject_v(N,E0,T,delta):
    
    xmin = 0
    xmax = E0
    n_accept = 0
    x_list = []
    while n_accept < N:
        t = (xmax - xmin) * np.random.rand() + xmin
        y = np.random.rand()
        if y < vaughan_pdf_corrected(t, E0, T,delta) / vaughan_pdf_corrected(T, E0, T,delta):
            n_accept += 1
            x_list.append(t)
    return x_list

def sey_coefficient(E, theta, E_th=0, Emax=500, delta_max=3, s=1, alpha=0.25):
    """
    Computes the SEY coefficient using the Modified Vaughan's model.
    
    Parameters:
    E (float): Impacting electron kinetic energy.
    theta (float): Incident angle with respect to the surface normal (in radians).
    E0 (float): Work function related parameter.
    Emax (float): Maximum energy for SEY.
    delta_max (float): Maximum SEY value.
    k_delta (float): Roughness factor related to delta.
    k_E (float): Roughness factor related to energy.
    alpha (float): Material-dependent fit parameter.
    delta_low (float): SEY value at low impacting energies (default is 1).

    theta = np.radians(0)  # Incident angle in radians
    E0 = 0  # Work function related parameter
    Emax = 340.0  # Maximum energy for SEY
    delta_max = 3.7  # Maximum SEY value
    k_delta = 1.0  # Roughness factor related to delta
    k_E = 1.0  # Roughness factor related to energy
    K_delta and K_E are roughness factors related to delta and energy and are now s
    alpha = 0.25  # Material-dependent fit parameter
    
    Returns:
    float: The SEY coefficient δ(E, θ).
    """
    Emax_theta = Emax * (1 + s * theta**2 / (2 * np.pi))
    
    delta_max_theta = delta_max * (1 + s * theta**2 / (2 * np.pi))
    
  
    v_E = (E - E_th) / (Emax_theta - E_th)
    if v_E < 0:
        delta = 0
    elif v_E < 1:
        delta = delta_max_theta * (v_E * np.exp(1 - v_E))**0.56
    elif 1 <= v_E <= 3.6:
        delta = delta_max_theta * (v_E * np.exp(1 - v_E))**alpha
    else:  # v(E) > 3.6
        delta = delta_max_theta * 1.125 * v_E**-0.35
    
    return delta




def sey_coefficient_guest(E, theta, Emax=500, delta_max=4, alpha=0.6, beta=0.65):
    C = np.cos(theta)*np.sqrt(E/ Emax)
    C = abs(C)
    #beta = 0.62 if E < Emax else 0.65
    term1 = (E / Emax) * np.sqrt(C)
    exp_term = np.exp(alpha * (1 - C) + beta * (1 - term1))
    return delta_max * term1 * exp_term

# Precompute the CDF and its inverse
def precompute_inverse_cdf(E0, T, delta, n_points=1000):
    E_values = np.linspace(0, E0, n_points)
    pdf_values = vaughan_pdf_corrected(E_values, E0, T, delta)
    cdf_values = np.cumsum(pdf_values) * (E_values[1] - E_values[0])  # Trapezoidal approximation
    cdf_values /= cdf_values[-1]  # Ensure CDF is normalized to 1
    
    # Ensure the CDF is strictly increasing
    cdf_values = np.clip(cdf_values, 0, 1)
    
    inverse_cdf_func = interp1d(cdf_values, E_values, bounds_error=False, fill_value=(0, E0))
    return inverse_cdf_func

# Generate random samples using the precomputed inverse CDF
def generate_samples_precomputed(n_samples, inverse_cdf_func):
    uniform_randoms = np.random.uniform(0, 1, n_samples)
    samples = inverse_cdf_func(uniform_randoms)
    return samples


def inverse_cdf_output(n_samples, E0, T,delta):
    inverse_cdf_func = precompute_inverse_cdf(E0, T,delta)
    samples = generate_samples_precomputed(n_samples, inverse_cdf_func)
    
    return samples




# Parameters
E0 = 150  # Example value for E0
T = 10  # Example value for T
n_samples = 100000
delta = 1
# Precompute the inverse CDF



plot = False

if plot == True:
    total_pos = []
    for i in range(0, 10):
        impact_energy = 600
        angle = 0
        pos_0 = test(impact_energy, angle)
        total_pos.append(pos_0)
    
    plt.hist(total_pos, bins=10, density=True, alpha=0.6)

    plt.figure()
    E0 = 150  # eV, assumed value
    T = 10  # eV, assumed temperature
    # Energy range for plotting
    E = np.linspace(0, 600, 5000)
    r = 0.025  # Radius of the cylinder
    d = 0.2  # Height of the cylinder
    end_position = [0.0025,     0.,         0.01779044]
    end_velocity =[0.86147291, 0.,          6.13108814]
    yield_curve = [sey_coefficient(E_i, theta=0, E_th=10, Emax=500, delta_max=4, k_delta=1, k_E=1, alpha=0.25) for E_i in E]
    
    theta_rad , theta = calculate_theta_cylinder(end_velocity, end_position, r)
    print(theta)
    yield_curve_rad = [sey_coefficient(E_i, theta=np.deg2rad(32), E_th=10, Emax=500, delta_max=4, k_delta=1, k_E=1, alpha=0.25) for E_i in E]
    plt.figure()
    plt.plot(E, yield_curve, 'b-', lw=2, label = f"theta = {0}")
    plt.plot(E, yield_curve_rad, 'r-', lw=2, label = f"theta = {32}")
    plt.title("SEY curve for a flat surface")
    plt.xlabel("Incident electron energy in eV")
    plt.ylabel("SEY coefficient")
    plt.grid(True)
    plt.legend()

    plt.figure()

    # Calculate the PDF values using the correct parameters
    pdf_values_corrected = vaughan_pdf_corrected(E, E0, T, delta)
    plt.plot(E, pdf_values_corrected, 'r-', lw=2)
    x_list= accept_reject_v(10000,E0,T,delta)
    samples = inverse_cdf_output(n_samples, E0, T,delta)
    # Plot the PDF
    plt.hist(x_list, bins=50, density=True, alpha=0.6, label='Accept Reject histogram')
    plt.hist(samples, bins=50, density=True, alpha=0.6, label='Precomputed histogram')


    # Plot the corrected PDF
    plt.plot(E, pdf_values_corrected, 'r-', lw=2)
    plt.title(f"PDF - Vaughan  E0 ={E0} eV, T= {T} eV")
    plt.xlabel("Secondary electron energy in eV")
    plt.ylabel("PDF - Vaughan")
    plt.grid(True)
    plt.legend()
    plt.show()


    sigmaFit = 0.0828
    muFit = 3.6621
    x = np.linspace(0.0001, 150, 1000)
    y = p(x,sigmaFit,muFit)
    N = 10000
    plt.plot(x,y, label = f"pdf")
    x_list = accept_reject(N,sigmaFit,muFit)
    bin_counts, bin_edges, patches = plt.hist(x_list, bins=100, density=True, alpha=0.6, label='Accept Reject histogram')

    #pdf_electrons  = pd.read_csv(f"pdf_python.txt", skiprows= 3, names = ["energy","dis"], sep="\t+")
    #plt.plot(pdf_electrons["energy"],pdf_electrons["dis"], label = f"pdf")
    plt.title("showing pdf v energy")
    plt.xlabel('energy(eV)')
    plt.ylabel('pdf')
    

   





#y = p(x,1.0828,1.66)
#plt.plot(x,y)

#plt.figure()
#y_g = np.linspace(0, 1, 1000)
#x_g= sc.gammaincinv(2, y_g)

#plt.plot(y_g,x_g)

#print(x_g)
#y_pdf = 2.92779*(x/(7.5**2))*np.exp(-x/7.5)
#plt.plot(x,y_pdf)



plt.show()