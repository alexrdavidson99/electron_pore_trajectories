
from matplotlib import pyplot as plt
import pandas as pd

voltages = [800, 1000, 1100, 1200,1300]

for i in range(len(voltages)):
    data = pd.read_csv(f"p_num_v_time_{voltages[i]}v.txt", skiprows= 5, names = ["time","number_of_electrons"], sep="\t+")
    plt.plot(data["time"],data["number_of_electrons"], label = f"{voltages[i]}V")


plt.title("showing number of electrons v time when changing voltage")
plt.xlabel('time(ns)')
plt.ylabel('number_of_electrons(n)')
plt.grid(True, which="both", ls="-")

plt.legend()
plt.show()
gain = pd.read_csv("pos_z_x_1300.txt", names=["x","y"], skiprows=3, sep="\t+")
print(len(gain["y"]))

