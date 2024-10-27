library(readr)
library(dplyr)
library(tidyr)
library(ggdendro)
library(VennDiagram)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(scales)

##########Venn diagram

# Function to format numbers with commas and percentages
format_with_commas_and_percent <- function(x, total) {
  formatted_number <- format(x, big.mark = ",", scientific = FALSE)
  percentage <- sprintf("(%.2f%%)", (x / total) * 100)
  return(paste(formatted_number, "\n", percentage))
}

# Define the sets and their intersections
mm10 <- 226339
D <- 8901818
B <- 3066205
mm10_D <- 8788344
mm10_B <- 14766109
D_B <- 9797161
mm10_D_B <- 15309954

# Calculate total elements
total_elements <- mm10 + D + B + mm10_D + mm10_B + D_B + mm10_D_B

# Calculate areas for each region
area1 <- mm10 + mm10_D + mm10_B + mm10_D_B
area2 <- D + mm10_D + D_B + mm10_D_B
area3 <- B + mm10_B + D_B + mm10_D_B
n12 <- mm10_D + mm10_D_B
n23 <- D_B + mm10_D_B
n13 <- mm10_B + mm10_D_B
n123 <- mm10_D_B

# Create the Venn diagram
venn.plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("mm10", "D", "B"),
  fill = c("#E6E6FA", "#FFB6C1", "#98FB98"),
  lty = "blank",
  cex = 0,
  cat.cex = 2,
  cat.col = c("purple", "pink", "green"),
  euler.d = TRUE,
  scaled = TRUE,
  alpha = 0.5,
  ind = TRUE
)

# Open a PNG graphics device
png("venn_diagram_bold_set_names_only.png", width = 2500, height = 2500, res = 300)

# Draw the base Venn diagram
grid.draw(venn.plot)

# Add custom labels with bold set names and regular numbers/percentages
add_label <- function(x, y, name, value, total) {
  grid.text(label = name, x = x, y = y + 0.03, gp = gpar(fontface = "bold", fontsize = 14))
  grid.text(label = format_with_commas_and_percent(value, total), x = x, y = y - 0.02, gp = gpar(fontsize = 12))
}

add_label(0.2, 0.8, "mm10", mm10, total_elements)
add_label(0.8, 0.8, "D", D, total_elements)
add_label(0.5, 0.2, "B", B, total_elements)
add_label(0.5, 0.75, "mm10&D", mm10_D, total_elements)
add_label(0.3, 0.5, "mm10&B", mm10_B, total_elements)
add_label(0.7, 0.5, "D&B", D_B, total_elements)
add_label(0.5, 0.5, "mm10&D&B", mm10_D_B, total_elements)
dev.off()

print("Venn diagram with bold set names only has been saved as 'venn_diagram_bold_set_names_only.png'")

##############################odgi similarity plot, not useful anymore
path_dist_tsv <- 'og_sim.tsv'

# Read sparse matrix and rename haplotypes
sparse_matrix_df <- read_tsv(path_dist_tsv) %>%
  mutate(across(c(group.a, group.b), ~ case_when(
    . == "REF" ~ "mm10",
    . == "DBA2J" ~ "D",
    . == "T2T_C57BL_6J" ~ "B",
    TRUE ~ .
  )))

# Print unique values to check renaming
print("Unique haplotype names after renaming:")
print(unique(c(sparse_matrix_df$group.a, sparse_matrix_df$group.b)))

# Prepare distance matrix
jaccard_dist_df <- sparse_matrix_df %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.similarity) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.similarity)

# Set row names and remove the group.a column
rownames(jaccard_dist_df) <- jaccard_dist_df$group.a
jaccard_dist_df <- jaccard_dist_df[, -1]

# Perform hierarchical clustering
jaccard_hc <- jaccard_dist_df %>% as.dist() %>% hclust()

# Create dendrogram object
dend <- as.dendrogram(jaccard_hc)
dend_data <- dendro_data(dend, type = "rectangle")

# Create the plot
p <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(dend_data), 
            aes(x = x, y = y, label = label, hjust = 0),
            size = 3, angle = 90) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_minimal() +
  labs(
    title = 'B-D-ref graph sequence similarity',
    x = 'Haplotype',
    y = 'Jaccard similarity'
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
ggsave("dendrogram_ggplot_renamed.png", p, width = 5, height = 9)


############################Nucleotide bases plot
# Deactivate scientific notation
options(scipen = 999)

# Create the data frame with N as the last entry
data <- data.frame(
  Base = factor(c("A", "C", "T", "G", "N"), levels = c("A", "C", "T", "G", "N")),
  Count = c(839712623, 590088595, 835116173, 591908088, 82077520)
)

# Calculate the percentage for each nucleotide base
data$Percentage <- (data$Count / sum(data$Count)) * 100

# Create the bar chart with percentages displayed on the bars
ggplot(data, aes(x = Base, y = Count, fill = Base)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            vjust = -0.5, size = 4) +  # Add percentage labels above the bars
  labs(title = "Nucleotide Base Counts for Sample D_C_mm10.gfa",
       x = "Nucleotide Base",
       y = "Count") +
  scale_y_continuous(labels = comma) +  # Format y-axis with commas
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

##############################################################Panel e, piechart plot

#setwd("/home/flavia/work/DBA2J/GWAS/debug_graph")
# Define the data
data <- data.frame(
  category = c("Nodes not duplicated", "Nodes duplicated", "Unmapped nodes"),
  count = c(14130807, 1277791, 9058896)
)

# Calculate percentages
data$percentage <- data$count / sum(data$count) * 100
# Define the data
data <- data.frame(
  category = c("Not duplicated", "Duplicated", "Unmapped"),
  count = c(14130807, 1277791, 9058896)
)

# Calculate percentages
data$percentage <- data$count / sum(data$count) * 100

# Create the plot
p <- ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "Not duplicated" = "#ae4371", 
    "Duplicated" = "#66c2a3",
    "Unmapped" = "#7570b3"
  )) +
  theme_void() +
  theme(legend.position = "bottom") +
  geom_text(aes(label = sprintf("%.1f%%\n(%s)", 
                                percentage, 
                                format(count, big.mark = ","))),
            position = position_stack(vjust = 0.5)) +
  labs(title = "Distribution of pangenome nodes",
       fill = "Category")

# Display the plot
ggsave("pangenomedistribution.png", p, width = 4, height = 4, dpi = 300)


############################Type of duplicated nodes
count_lines <- function(file_path) {
  length(readLines(file_path))
}
diff_pos_count <- count_lines("unique_diff_pos.txt")
diff_chr_count <- count_lines("unique_diff_chr.txt")
shared_count <- count_lines("shared_nodes_diff_chr.txt")

# Read length data
length_data <- read.table("uniqduplicates_with_lengths.txt", col.names = c("node", "length"))

# Create data frame
data_duplicates <- data.frame(
  category = c("different pos", "different chrm", "both"),
  count = c(1231772, 32689, 13330)
)

# Calculate total count
total_count <- sum(data_duplicates$count)

# Set factor levels for ordering
data_duplicates$category <- factor(data_duplicates$category, 
                                   levels = c("both", "different chrm", "different pos"))

# Calculate percentages correctly
data_duplicates <- data_duplicates %>%
  mutate(
    percentage = (count / total_count) * 100
  )



options(scipen = 999)

# Create the plot for duplicates (count)
p2 <- ggplot(data_duplicates, aes(x = " ", y = count, fill = category)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf),
            fill = "#66c2a3", alpha = 0.2) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +  # Cambiato da 0.7 a 0.4
  scale_fill_manual(values = c("different pos" = "#4372ae", "different chrm" = "#ae4372", 
                               "both" = "#a7ae43")) +
  labs(title = " ",
       x = "Duplicated nodes",
       y = "Count",
       fill = " ") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) + 
  scale_y_continuous(labels = scales::comma) +coord_flip()
# Ensure the category is a factor with the desired order
data_duplicates$category <- factor(data_duplicates$category, 
                                   levels = c("both", "different chrm", "different pos"))

# Calculate percentages correctly
total_count <- sum(data_duplicates$count)
data_duplicates <- data_duplicates %>%
  mutate(
    percentage = count / total_count * 100
  )

options(scipen = 999)

# Create the plot for duplicates (percentage)
p3 <- ggplot(data_duplicates, aes(x = " ", y = percentage, fill = category)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100),
            fill = "#66c2a3", alpha = 0.2) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +
  scale_fill_manual(values = c("both" = "#a7ae43", 
                               "different chrm" = "#ae4372", 
                               "different pos" = "#4372ae"),
                    breaks = c("different pos", "different chrm", "both")) +
  labs(title = " ",
       x = "Duplicated nodes",
       y = "Percentage",
       fill = " ") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.y = element_blank(),  # Changed from x to y due to coord_flip()
    axis.ticks.y = element_blank(),  # Changed from x to y due to coord_flip()
    panel.grid.major.y = element_blank(),  # Changed from x to y due to coord_flip()
    panel.grid.minor.y = element_blank(),  # Changed from x to y due to coord_flip()
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     breaks = seq(0, 100, by = 20),
                     expand = expansion(mult = c(0, 0.05))) +  # Added expand for better spacing
  coord_flip()

print(p2)

################ The histogram plot of duplicated nodes
p4 <- ggplot(length_data, aes(x = length)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#66c2a0", alpha = 0.1)+
  geom_histogram(binwidth = 50) +
  labs(title = " ",
       x = "Duplicated node length",
       y = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_y_log10() + xlim(0,10000) 


###################################Lenght nodes by different type plot
pos_data <- read.table("pos_diff_with_length.txt", header = FALSE, col.names = c("Node", "Length"))
chr_data <- read.table("chr_diff_with_length.txt", header = FALSE, col.names = c("Node", "Length"))
shared_data <- read.table("shared_with_length.txt", header = FALSE, col.names = c("Node", "Length"))

# Add a category column
pos_data$Category <- "different pos"
chr_data$Category <- "different chrm"
shared_data$Category <- "both"

# Combine the data into one data frame
combined_data <- bind_rows(pos_data, chr_data, shared_data)
combined_data$Category <- factor(combined_data$Category, 
                                   levels = c("both","different chrm", "different pos"))
# Create a new column for length categories
combined_data$Length_Category <- cut(
  combined_data$Length,
  breaks = c(-Inf, 100, 1000, 10000, Inf),
  labels = c("0-100", "100-1000", "1000-10000", ">10000")
)

# Calculate the total count for each Length_Category
total_counts <- combined_data %>%
  group_by(Length_Category) %>%
  summarise(Total_Length_Category = n(), .groups = 'drop')

# Summarize the data by length category and category type, then join with the total counts
summary_data <- combined_data %>%
  group_by(Length_Category, Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  left_join(total_counts, by = "Length_Category") %>%
  mutate(Percentage = Count / Total_Length_Category)  # Calculate the percentage per group

# Plot p5 with percentages ranging from 0 to 100 on the y-axis
p5 <- ggplot(summary_data, aes(x = Length_Category, y = Count, fill = Category)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "#66c2a0", alpha = 0.1) +
  theme(legend.position = "bottom") +
  geom_bar(stat = "identity", position = "stack", width = 0.3) +
  scale_fill_manual(values = c("both" = "#a7ae43", 
                               "different chrm" = "#ae4372", 
                               "different pos" = "#4372ae"),
                    breaks = c("different pos", "different chrm", "both")) +
  labs(title = " ",
       x = "Duplicated node length",
       y = "Count",
       fill = " ") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(5.5, 20, 5.5, 5.5)) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  coord_flip()



# Ensure the Percentage is calculated as a value between 0 and 100
summary_data <- summary_data %>%
  mutate(Percentage = Percentage * 100)

# Plot p5 with percentages ranging from 0 to 100 on the y-axis
p5bis <- ggplot(summary_data, aes(x = Length_Category, y = Percentage, fill = Category)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100),
            fill = "#66c2a3", alpha = 0.1) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("different pos" = "#4372ae", 
                               "different chrm" = "#ae4372", 
                               "shared" = "#a7ae43")) +
  labs(title = "Duplicated nodes categorized by length",
       x = "Length Nodes",
       y = "Percentage",
       fill = "Category") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%"))

##########plot overlap/non overlap repeats
overlap <- read_tsv("overlap_dup_with_rmsk.txt", col_names = c("chr", "start", "end", "repeat_class"))
non_overlap <- read_tsv("non_overlap_dup_regions.txt", col_names = c("chr", "start", "end"))

# Calculate total counts
total_count <- nrow(overlap) + nrow(non_overlap)

# Calculate percentages
overlap_percent <- nrow(overlap) / total_count * 100
non_overlap_percent <- nrow(non_overlap) / total_count * 100

# Create a summary dataframe
summary_data <- data.frame(
  Category = c("Overlap", "Non-overlap"),
  Percentage = c(overlap_percent, non_overlap_percent)
)

# First, ensure Category is a factor and set the levels in the desired order
summary_data$Category <- factor(summary_data$Category, levels = c("Overlap", "Non-overlap"))


# Group by major families
repeat_summary <- overlap %>%
  mutate(major_family = sapply(strsplit(as.character(repeat_class), "/"), `[`, 1)) %>%
  group_by(major_family) %>%
  summarise(count = n()) %>%
  mutate(Percentage = (count / sum(count)) * 100,
         Category = "Overlap") %>%
  arrange(desc(count))  # Sort by frequency

# First, let's verify the exact data in summary_data
print("Categories in summary_data:")
print(unique(summary_data$Category))

# Create the basic plot first to see what categories are being used
basic_plot <- ggplot(summary_data, aes(x = "", y = Percentage, fill = Category))
print("Categories being used in plot:")
print(levels(basic_plot$data$Category))

# Create color mappings
blue_shades <- c(
  "#00008B",  # Dark blue (for LINE)
           "#0066CC",  # Medium-dark blue
           "#1E90FF",  # Dodger blue
           "#00BFFF",  # Deep sky blue
           "#87CEEB",  # Sky blue
           "#B0E0E6",  # Powder blue
           "#ADD8E6",  # Light blue
           "#4169E1",  # Royal blue
           "#6495ED",  # Cornflower blue
           "#A1CAF1",  # Baby blue
           "#89CFF0",  # Carolina blue
           "#F0F8FF",  # Alice blue
           "#B9D9EB",  # Columbia blue
           "#CBDAF2"   # Beau blue
)

# Let's try forcing the Category to be a factor with specific levels
summary_data$Category <- factor(summary_data$Category, 
                                levels = c("Non-overlap", "Overlap"))

# Create the colors with explicit mapping
repeat_colors <- setNames(blue_shades[1:nrow(repeat_summary)], repeat_summary$major_family)
color_mapping <- c("Non-overlap" = "#2BECCC")  # Start with Non-overlap
color_mapping <- c(color_mapping, repeat_colors)  # Add repeat colors

# Create the plot with forced mappings
# Create the plot with forced mappings
p6 <- ggplot(summary_data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_bar(data = repeat_summary, 
           aes(x = "", y = Percentage * (summary_data$Percentage[summary_data$Category == "Overlap"]/100), 
               fill = major_family),
           stat = "identity", width = 0.4) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100), 
            fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 100),
            fill = "#66c2a3", alpha = 0.2) +
  scale_fill_manual(
    values = color_mapping,
    limits = c("Non-overlap", repeat_summary$major_family),
    guide = guide_legend(nrow = 2)  # Force the legend to be on one row
  ) +
  labs(title = " ",
       x = "Duplicated nodes",
       y = "Percentage",
       fill = "Repeats") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.spacing.x = unit(0.3, 'cm'),  # Adjusts spacing between legend items
    legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), 
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 100, by = 20)) +
  coord_flip()

# Plot the figure
print(p6)


# Print debugging information
print("Final color mapping:")
print(color_mapping)
print("Category levels in summary_data:")
print(levels(summary_data$Category))

# Combine the plots
combined_plot <- (p3 | p6 ) / (p5 | p4) + plot_layout(heights = c(1, 1,1))

# Display the combined plot
print(combined_plot)
ggsave("combinedplot.png", combined_plot, width = 14, height = 7, dpi = 300)
