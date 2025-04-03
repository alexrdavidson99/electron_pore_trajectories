import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_data(file_name):
    return pd.read_csv(fr"./data_from_ion_pot/{file_name}",names=['time'])
   
def plot_data_outputs(data):
    n, bins = np.histogram(data['time'], bins=100)  # density=True)
    n_max = n.max()
    n_normalized = n / n_max  # Normalize the bin counts
    return n_normalized, bins

list_of_files = ['times_by_energy_voltage_uni.csv', 'times_by_energy_voltage_ar.csv', 'times_by_energy_voltage_cdf.csv']
list_of_files = ['times_by_energy_voltage_cdf.csv']
for file_name in list_of_files:
    data = read_data(file_name)
    n_normalized, bins = plot_data_outputs(data)
    split_name = file_name.split('_')
    plt.bar(bins[:-1], n_normalized, width=(bins[1] - bins[0]), alpha=0.5, label=split_name[4])

plt.xlabel('time (s)')
plt.ylabel('normalized count')
plt.title('time taken for ions to reach the pore')
plt.legend()
plt.show()




