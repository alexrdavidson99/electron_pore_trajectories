import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mticker
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad

# My axis should display 10⁻¹ but you can switch to e-notation 1.00e+01
def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"  # remove int() if you don't use MaxNLocator
    # return f"{10**val:.2e}"      # e-Notation

def model_func(x, a, b):
    return a * np.exp(b * x) 


def accept_reject(N,a,b,xmin,xmax):
    
    
    n_accept = 0
    x_list = []
    
    while n_accept < N:
        t = (xmax - xmin) * np.random.rand() + xmin
        y = np.random.rand()
        if y < model_func(t,a,b) / model_func(xmax,a,b):
            n_accept += 1
            x_list.append(t)
    return x_list

def normalize_theoretical_distribution(a, b, xmin, xmax):
    # Integrate the model function over the range [xmin, xmax]
    integral, _ = quad(lambda x: model_func(x, a, b), xmin, xmax)
    
    # Normalize the model function by dividing by the integral
    def normalized_model_func(x):
        return model_func(x, a, b) / integral
    
    return normalized_model_func


def sample_start_position(cdf, positions):
    random_value = np.random.rand()
    index = np.searchsorted(cdf, random_value)
    return positions[index]

# Read the CSV files
df_1ev = pd.read_csv('prob_results_1ev.csv')
df_2ev = pd.read_csv('prob_results_2ev.csv')
df_3ev = pd.read_csv('prob_results_3ev.csv')

# Extract distance values and count results
distance_values = 46e-2-df_1ev['distance_values']
count_results_1ev = np.log10(df_1ev.drop(columns=['distance_values']).values)
count_results_2ev = np.log10(df_2ev.drop(columns=['distance_values']).values)
count_results_3ev = np.log10(df_3ev.drop(columns=['distance_values']).values)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot data for 1 eV
for i, voltage in enumerate(df_1ev.columns[1:]):
    ax.plot(distance_values, count_results_1ev[:, i], zs=1, zdir='z', label=f'1 eV, {voltage}')

# Plot data for 2 eV
for i, voltage in enumerate(df_2ev.columns[1:]):
    ax.plot(distance_values, count_results_2ev[:, i], zs=2, zdir='z', label=f'2 eV, {voltage}')

# Plot data for 3 eV
for i, voltage in enumerate(df_3ev.columns[1:]):
    ax.plot(distance_values, count_results_3ev[:, i], zs=3, zdir='z', label=f'3 eV, {voltage}')
    
# Set labels
ax.set_xlabel('Distance Values (mm)')
ax.set_ylabel('Probability of exiting the pore')

ax.set_zlabel('Energy (eV)')
ax.set_title('Ion escape probability for Different Energies and Voltages')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax.legend()



#reading in end position data
end_position = pd.read_csv('end_positions_2.txt', names=["z","hits"], sep=",")
end_position_norm = end_position["hits"]/end_position["hits"].max()

popt, _ = curve_fit(model_func, end_position["z"], end_position_norm)
a = popt[0]
b = popt[1]
xmin = 0.04
xmax = end_position["z"].max()
N = 10000
print(f"y = {a} * exp({b} * x)")
data_from_accept_reject = accept_reject(N,a,b,xmin,xmax)
fit_x = np.linspace(xmin, xmax, 1000)
fit_y = model_func(fit_x, *popt)

normalized_model_func = normalize_theoretical_distribution(a, b, xmin, xmax)
normalized_model = np.array([normalized_model_func(x) for x in fit_x])
sampled_positions_df = pd.DataFrame(data_from_accept_reject, columns=["start_position"])
sampled_positions_df.to_csv("sampled_start_positions_ar.csv", index=False)


plt.figure()
plt.scatter(end_position["z"], end_position_norm, alpha=0.5)
fit_y = model_func(fit_x, *popt) / np.max(model_func(fit_x, *popt))
plt.plot(fit_x, fit_y, 'r-', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
plt.plot(fit_x, normalized_model, label='Normalized Theoretical Distribution', color='r', linewidth=2)
plt.hist(data_from_accept_reject, bins=int(N/25), alpha=0.75,density=True, label='Accept Reject histogram')

#plt.yscale('log')
plt.title("Impact position_z_1000_runs")
plt.xlabel("position (mm)")
plt.ylabel("counts_norm")
plt.legend()
plt.figure()




#matched_indices = np.searchsorted(distance_values, end_position["z"])


one_ev_not_loged = df_1ev.drop(columns=['distance_values']).values
plt.scatter(distance_values, one_ev_not_loged[:,2], alpha=0.5,label="1ev ion escape probability")
plt.scatter(end_position["z"], end_position_norm, alpha=0.5, label= "normalised number of electrons hitting the pore wall")
plt.ylabel("probability")
plt.xlabel("distance from the pore entrance (mm)")
plt.legend()
print(one_ev_not_loged[:,1][:-1])
print(end_position_norm[::-1])

plt.figure()

for i, voltage in enumerate(df_1ev.columns[1:]):
    matched_count_results_1ev = one_ev_not_loged[:,i][:-1] * end_position_norm[::-1]
    plt.scatter(end_position["z"][::-1], matched_count_results_1ev, alpha=0.5, label=f"total 1ev ion escape probability at {voltage} V")
    print(f'V{voltage}')
    if voltage == "V_500":
        plotted_data_df = pd.DataFrame({
            "z": end_position["z"][::-1],
            "matched_count_results_1ev": matched_count_results_1ev
        })
    plotted_data_df.to_csv("plotted_data.csv", index=False)


plt.ylabel("probability")
plt.xlabel("distance from the pore entrance (mm)")
plt.yscale('log')
plt.legend()

probability_distribution = matched_count_results_1ev / np.sum(matched_count_results_1ev)
cdf = np.cumsum(probability_distribution)


sampled_positions = [sample_start_position(cdf, end_position["z"][::-1]) for _ in range(10000)]
#sampled_positions = [pos for pos in sampled_positions if -0.4575 <= pos <= -0.0328]
# Plot the sampled start positions
plt.figure()
sampled_positions = np.array(sampled_positions)
plt.hist(sampled_positions, bins=50, alpha=0.75, label='Sampled Start Positions')
plt.xlabel("Start Position (mm)")
plt.ylabel("Frequency")
plt.yscale('log')
plt.xlim(0,0.5)
plt.title("Histogram of Sampled Start Positions from CDF")
plt.legend()

sampled_positions_df = pd.DataFrame(sampled_positions, columns=["start_position"])
sampled_positions_df.to_csv("sampled_start_positions_CDF_likelyhood.csv", index=False)


# Show plot
plt.show()