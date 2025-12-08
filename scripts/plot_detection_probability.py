
import numpy as np
import matplotlib.pyplot as plt
import os

#
# Parameters
#

# NOTE: Write the prefix-name of the input files you want to plot here
configs = [
    ("single_slit_experiment",  "Single slit"),
    ("double_slit_experiment",  "Double slit"),
    ("triple_slit_experiment",  "Triple slit"),
]

M  = 201             # total grid points per dimension
dt = 2.5e-5          # delta time
T  = 0.002           # Time interval
Nt = int(T / dt) + 1 # should be 81

# Detection probability at x = 0.8 at t = 0.002
x_screen = 0.8
t_target = 0.002

h = 1.0 / (M - 1)                   # consistent with C++ grid
i_screen = int(round(x_screen / h)) # should be 160 for M=201
n_t = int(round(t_target / dt))     # should be 80 for dt=2.5e-5

#
# Get path to the input files
#

folder_path = "output/"
if not os.path.isdir(folder_path):
    folder_path = "../" + folder_path

while not os.path.isdir(folder_path):
    folder_path = input("Please enter the path to the folder containing the input files: ")

#
# Plot
#

y = np.linspace(0.0, 1.0, M)

for prefix, label in configs:
    # Load Data
    file_name = os.path.join(folder_path, prefix + "_wave.bin")
    data = np.fromfile(file_name, dtype=np.complex128)
    U = data.reshape((M, M, Nt), order="F")

    # Extract the 1D wave function along y at the screen
    u_line = U[i_screen, :, n_t]
    # Unnormalised probabilities along the screen
    q = np.abs(u_line)**2
    # Normalized detection probability distribution along the screen
    # (conditional on detection at x = i_screen)
    p_cond = q / q.sum()

    plt.plot(y, p_cond, label=label)

plt.xlabel(r"$y$")
plt.ylabel(r"$p(y\,|\,x=0.8;\,t=0.002)$")
plt.title("Detection probability at screen x = 0.8")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(folder_path, "detprob_single_double_triple.pdf"))
plt.close()
