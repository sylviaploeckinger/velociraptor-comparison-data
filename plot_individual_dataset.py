"""
Plots an individual dataset and shows the result.
"""

import sys
import matplotlib.pyplot as plt
from velociraptor.observations import load_observation

obs = load_observation(sys.argv[1])

fig, ax = plt.subplots()
ax.loglog()

obs.plot_on_axes(ax)

ax.set_xlabel(f"{obs.x_description} $\\left[{obs.x_units.latex_repr}\\right]$")
ax.set_ylabel(f"{obs.y_description} $\\left[{obs.y_units.latex_repr}\\right]$")
ax.set_title(obs.name)

ax.legend()
fig.tight_layout()

ax.text(
    0.025,
    0.025,
    obs.comment,
    ha="left",
    va="bottom",
    transform=ax.transAxes,
    wrap=True,
    fontsize=8,
)

plt.show()
