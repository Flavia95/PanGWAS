import matplotlib.pyplot as plt

# Data for the bar chart
bases = ['A', 'C', 'T', 'G', 'N']
counts = [800_000_000, 560_000_000, 795_000_000, 560_000_000, 80_000_000]
percentages_bar = [28.6, 20.1, 28.4, 20.1, 2.8]

# Set up figure with 1:1.7 aspect ratio
fig, ax = plt.subplots(figsize=(8.5, 5))  # 8.5:5 is close to 1:1.7 aspect ratio

# Create the bar plot
bars = ax.bar(bases, [count / 1_000_000 for count in counts], color=['salmon', 'olive', 'mediumseagreen', 'cornflowerblue', 'orchid'])

# Annotate with percentages
for bar, percentage in zip(bars, percentages_bar):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f'{percentage}%',
        ha='center',
        va='bottom',
        fontsize=12,
        fontweight='bold'
    )

# Labeling for the bar chart
ax.set_xlabel("Nucleotide Base", fontsize=18)
ax.set_ylabel("Counts (M)", fontsize=18)
ax.tick_params(axis='both', labelsize=16)  # Increase font size for tick labels
ax.yaxis.set_ticks(range(0, int(max(counts) / 1_000_000) + 400, 400))  # Set y-ticks every 400M

# Save the figure with transparent background
plt.savefig("barplot_vertical_horizontal.png", format="png", dpi=300, transparent=True)
plt.savefig("barplot_vertical_horizontal.svg", format="svg", dpi=300, transparent=True)

# Display the plot
plt.show()
