import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import pandas as pd
import random


def p(x,sigmaFit,muFit):
    # return x/(5**2)*np.exp(-(x**2)/(2*(5**2)))
    #sigmaFit = 1.0828
    #muFit = 1.6636
    #return 2.92779 * (x / (7.5 ** 2)) * np.exp(-x / 7.5)
    return (1 / (x * sigmaFit * (2 * np.pi) ** (1 / 2))) * np.exp(-(((np.log(x) - muFit) ** 2) / (2 * (sigmaFit ** 2))))


def accept_reject(N):
    xmin = 0
    xmax = 150
    #pmax = 0.12130
    pmax = 0.12544838348768622

    n_accept = 0
    x_list = []
    while n_accept < N:
        t = (xmax - xmin) * np.random.rand() + xmin
        y = np.random.rand()
        if y < p(t,1.0828,1.6636) / pmax:
            n_accept += 1
            x_list.append(t)
    return x_list


#bin_counts, bin_edges, patches = plt.hist(x, bins=100, density=True, alpha=0.6, label='Accept Reject histogram')

#x = np.linspace(0, 150, 1000)

#y = p(x,1.0828,1.6636)
#plt.plot(x,y, label = f"pdf")

#pdf_electrons  = pd.read_csv(f"pdf_python.txt", skiprows= 3, names = ["energy","dis"], sep="\t+")
#plt.plot(pdf_electrons["energy"],pdf_electrons["dis"], label = f"pdf")
#plt.title("showing pdf v energy")
#plt.xlabel('energy(eV)')
#plt.ylabel('pdf')
#plt.grid(True, which="both", ls="-")
#plt.legend()
#bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

# Generate some dummy error values, of the same dimensions as the bin counts
#y_error = np.random.rand(bin_counts.size)*0.005

# Plot the error bars, centred on (bin_centre, bin_count), with length y_error
#plt.errorbar(x=bin_centres, y=bin_counts,
#            yerr=y_error, fmt='o', capsize=2)





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