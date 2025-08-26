import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

# Load and preprocess data
file_path = 'simulation_results/parameter_estimates.csv'
data = pd.read_csv(file_path)
data["mean"] = pd.to_numeric(data["mean"])
data["SD"] = pd.to_numeric(data["SD"])

print(len(data))
# Filter out unrealistic or erroneous values
data = data[(data['mean'] > 0) & (data['mean'] <= 100000)]
data = data[(data['SD'] > 0) & (data['SD'] <= 100000)]

# Linear regression
slope, intercept, r_value, p_value, std_err = linregress(data['mean'], data['SD'])

# Font and plot settings
plt.rcParams.update({
    'font.family': 'Sans-Serif',
    'font.size': 12,
    'axes.linewidth': 1,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.direction': 'out',
    'ytick.direction': 'out'
})

# Create the plot
fig, ax = plt.subplots(figsize=(7, 4.5))

# Scatter and regression line
ax.scatter(data['mean'], data['SD'], color='#56B4E9', alpha=0.6, s=30)
ax.plot(data['mean'], intercept + slope * data['mean'], color='red', linewidth=2)

# Annotate regression equation
eq_text = f'y = {slope:.2f}x + {intercept:.2f}'
ax.text(0.05, 0.90, eq_text, transform=ax.transAxes, fontsize=12, color='red')

# Updated axis labels using LaTeX-style math
ax.set_xlabel(r'$\mu$ (bp)', labelpad=10)
ax.set_ylabel(r'$\sigma$ (bp)', labelpad=10)

# Limits and ticks
ax.set_xlim([0, 50000])
ax.set_ylim([0, 50000])
ax.set_xticks([0, 10000, 20000, 30000, 40000, 50000])
ax.set_yticks([0, 10000, 20000, 30000, 40000, 50000])

# Clean visual
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', length=6, width=1)

# Final layout and save
plt.tight_layout()
plt.savefig("simulation_results/mean_vs_sd.png", dpi=600)
plt.savefig("simulation_results/mean_vs_sd.pdf")

plt.close()