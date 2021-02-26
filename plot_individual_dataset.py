"""
Plots an individual dataset and shows the result.
"""

import sys
import matplotlib.pyplot as plt
import unyt
from velociraptor.observations import load_observations

try:
    lower = float(sys.argv[2])
    upper = float(sys.argv[3])
except:
    print("Redshift range not found, assuming 0.0 -> 1000.0")
    lower = 0.0
    upper = 1000.0

obs = load_observations(sys.argv[1], [lower, upper])

with unyt.matplotlib_support:
    fig, ax = plt.subplots()
    ax.loglog()

    for observation in obs:
        observation.plot_on_axes(ax)
        ax.set_title(observation.name)

ax.legend()
fig.tight_layout()

ax.text(
    0.025,
    0.025,
    obs[0].comment,
    ha="left",
    va="bottom",
    transform=ax.transAxes,
    wrap=True,
    fontsize=8,
)

plt.show()
