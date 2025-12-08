
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

#
# Simulation / file parameters
#

# NOTE: Change this prefix to plot single and tripple slit experiments
file_name_prefix = "double_slit_experiment"

M  = 201             # total grid points per dimension (must match main.cpp)
dt = 2.5e-5          # must match main.cpp
T  = 0.002           # total simulated time (adjust if needed)
Nt = int(T / dt) + 1 # number of time steps (must match main.cpp)

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
# Set up grid and plotting ranges
#

x = np.linspace(0.0, 1.0, M)
y = np.linspace(0.0, 1.0, M)
extent = [0.0, 1.0, 0.0, 1.0]  # [xmin, xmax, ymin, ymax]

# Time array (for convenience / labels)
times = np.arange(Nt) * dt

#
# Animation setup (inspired by the following URL)
# URL: https://anderkve.github.io/FYS3150/book/projects/project5.html#code-snippets
#

fontsize = 12

# Use probability density as z(x,y,t) = |u|^2, normalized per frame
def prob_frame(idx):
    u = U[:, :, idx]
    p = np.abs(u)**2
    p /= p.max() # normalize per frame so color scale is always [0, 1]
    return p

# Initial frame
z0 = prob_frame(0)

fig = plt.figure()
ax  = plt.gca()

# Color normalization based on initial frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z0))

# imshow wants shape (Ny, Nx), so we transpose and set origin='lower'
img = ax.imshow(
    z0.T,
    extent=extent,
    origin="lower",
    aspect="equal",
    cmap="viridis",
    norm=norm
)

ax.set_xlabel("$x$", fontsize=fontsize)
ax.set_ylabel("$y$", fontsize=fontsize)
ax.set_title("Probability density $|u(x,y,t)|^2$", fontsize=fontsize)

cbar = fig.colorbar(img, ax=ax)
cbar.set_label("Normalised $|u|^2$", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Text overlay showing current time
time_txt = ax.text(
    0.95, 0.95,
    f"t = {times[0]:.3e}",
    color="white",
    horizontalalignment="right",
    verticalalignment="top",
    transform=ax.transAxes,
    fontsize=fontsize
)

#
# Animation function: updates each frame
#

# You can skip frames to speed up the animation, e.g. frame_step = 2 or 5
frame_step = 2
frame_indices = np.arange(0, Nt, frame_step)

def update(frame_idx):
    i = frame_indices[frame_idx]

    z = prob_frame(i)
    # Option 1: keep fixed norm [0, 1] (already normalized)
    # Option 2: recompute norm per frame; here we keep [0, 1] so just set_data
    img.set_data(z.T)

    # Update time label
    time_txt.set_text(f"t = {times[i]:.3e}")

    return img, time_txt

anim = FuncAnimation(
    fig,
    update,
    frames=len(frame_indices),
    interval=30, # milliseconds between frames (adjust playback speed)
    blit=False
)

plt.tight_layout()
plt.show()

# Save the animation as an mp4
# (requires ffmpeg installed on system).
outfile = os.path.join(folder_path, file_name_prefix + "_animation.mp4")
anim.save(outfile, writer="ffmpeg", bitrate=10000, fps=30)
