
import numpy as np
import matplotlib.pyplot as plt
import os

#
# Parameters
#

# NOTE: Change this prefix to plot single and tripple slit experiments
file_name_prefix = "double_slit_experiment"

M  = 201             # total grid points per dimension
dt = 2.5e-5          # Delta time
T  = 0.002           # Time interval
Nt = int(T / dt) + 1 # should be 81

# Global Font Size
plt.rcParams.update({
    "font.size": 16,        # base font size
    "axes.titlesize": 18,   # title size
    "axes.labelsize": 16,   # x/y label size
    "legend.fontsize": 14,  # legend text
    "xtick.labelsize": 14,
    "ytick.labelsize": 14
})

#
# Get path to the input files
#

folder_path = "output/"
if not os.path.isdir(folder_path):
    folder_path = "../" + folder_path

while not os.path.isdir(folder_path):
    folder_path = input("Please enter the path to the folder containing the input files: ")

file_name = os.path.join(folder_path, file_name_prefix + "_wave.bin")

#
# Load data
#

# load cube: shape (M, M, Nt)
data = np.fromfile(file_name, dtype=np.complex128)
U = data.reshape((M, M, Nt), order="F")

#
# Plots
#

# time indices for t = 0, 0.001, 0.002
times   = [0.0, 0.001, 0.002]
indices = [int(t / dt) for t in times] # [0, 40, 80]

x = np.linspace(0.0, 1.0, M)
y = np.linspace(0.0, 1.0, M)
extent = [0.0, 1.0, 0.0, 1.0] # [xmin, xmax, ymin, ymax]

def plot_colormap(z, title, cbar_label, filename):
    plt.figure(figsize=(5, 4))
    im = plt.imshow(
        z.T,         # transpose so x is horizontal, y vertical
        origin="lower", 
        extent=extent,
        aspect="equal",
        cmap="viridis"
    )
    plt.xlabel("$x$-coordinate")
    plt.ylabel("$y$-coordinate")
    plt.title(title)
    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

for t, idx in zip(times, indices):
    u = U[:, :, idx]

    #
    # probability density p(x,y;t) = |u|^2
    #
    
    p = np.abs(u)**2
    
    # Rescale colour scale per time step
    use_direct_scaling = True # Toggle between direct normalization and sqrt scaling
    if use_direct_scaling:
        plot_colormap(
            p / p.max(), # normalised to max=1 at each t
            f"$t={t:.3e}$",
            "Normalized $|u|^2$",
            os.path.join(folder_path, f"prob_t{t:.3f}.pdf")
        )
    else:
        plot_colormap(
            np.sqrt(p / p.max()), # normalized with sqrt to max=1 at each t
            f"$t={t:.3e}$",
            "Normalized $\sqrt{|u|^2}$",
            os.path.join(folder_path, f"prob_t{t:.3f}.pdf")
        )

    #
    # Real part
    #
    
    plot_colormap(
        u.real,
        f"$t={t:.3e}$",
        "Re$(u)$",
        os.path.join(folder_path, f"re_u_t{t:.3f}.pdf")
    )

    #
    # Imaginary part
    #
    
    plot_colormap(
        u.imag,
        f"$t={t:.3e}$",
        "Im$(u)$",
        os.path.join(folder_path, f"im_u_t{t:.3f}.pdf")
    )
