import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Define the values for each section
n1 = 226339        # mm10 only
n2 = 8901818       # D only
n3 = 3066205       # B only
n12 = 8788344      # mm10_D (intersection of mm10 and D only)
n13 = 14766109     # mm10_B (intersection of mm10 and B only)
n23 = 9797161      # D_B (intersection of D and B only)
n123 = 15309954    # mm10_D_B (intersection of mm10, D, and B)

# Calculate the total to convert values to percentages
total = n1 + n2 + n3 + n12 + n13 + n23 + n123

# Calculate percentages
percentages = {
    '100': (n1 / total) * 100,
    '010': (n2 / total) * 100,
    '001': (n3 / total) * 100,
    '110': (n12 / total) * 100,
    '101': (n13 / total) * 100,
    '011': (n23 / total) * 100,
    '111': (n123 / total) * 100
}

# Actual values dictionary
actual_values = {
    '100': n1,
    '010': n2,
    '001': n3,
    '110': n12,
    '101': n13,
    '011': n23,
    '111': n123
}

# Create the Venn diagram with both actual numbers and percentages
plt.figure(figsize=(10, 8), dpi=300)  # Increase DPI for high resolution
venn = venn3(subsets=(percentages['100'], percentages['010'], percentages['110'], percentages['001'], percentages['101'], percentages['011'], percentages['111']), set_labels=('mm10', 'D', 'B'))

# Update colors
venn.get_patch_by_id('100').set_color("skyblue")
venn.get_patch_by_id('010').set_color("pink")
venn.get_patch_by_id('001').set_color("mediumorchid")

# Update labels with actual values and percentages
for subset, percent in percentages.items():
    if venn.get_label_by_id(subset) is not None:
        actual = actual_values[subset]
        venn.get_label_by_id(subset).set_text(f"{actual:,}\n({percent:.1f}%)")
        venn.get_label_by_id(subset).set_fontsize(16)  # Increased font size

# Set category labels in bold and larger font size
for label in venn.set_labels:
    label.set_fontsize(20)
    label.set_fontweight('bold')

# Save the image with a transparent background
plt.savefig("venn_diagram.png", format="png", dpi=300, transparent=True)
plt.savefig("venn_diagram.svg", format="svg", dpi=300, transparent=True)

# Display the plot without title
plt.show()
