
import os
import numpy as np
import matplotlib.pyplot as plt

#
# Get path to the input files
#

file_name_free_prop = "free_propagation_barrier_off_probability.csv"
file_name_double_slit = "double_slit_barrier_on_probability.csv"

folder_path = "output/"
if not os.path.isdir(folder_path):
    folder_path = "../" + folder_path

while not os.path.isdir(folder_path):
    folder_path = input("Please enter the path to the folder containing the input files: ")

#
# Load data
#

t_free, dev_free   = np.loadtxt(folder_path + file_name_free_prop, delimiter=",", unpack=True)
t_dslit, dev_dslit = np.loadtxt(folder_path + file_name_double_slit, delimiter=",", unpack=True)

#
# Plot absolute deviation on a log scale
#

plt.figure()
plt.semilogy(t_free,  np.abs(dev_free),   label="Free propagation ($v_0 = 0$)")
plt.semilogy(t_dslit, np.abs(dev_dslit),  label="Double-slit barrier ($v_0 = 10^{10}$)")

plt.xlabel("Time $t$")
plt.ylabel("Deviation from unity $|P(t) - 1|$")
plt.legend(title="Simulation cases", loc="best", frameon=True)
plt.grid(True, which="both", ls=":")
plt.tight_layout()

output_name = os.path.join(folder_path, "probability_deviation_plot.pdf")
plt.savefig(output_name)
plt.show()
