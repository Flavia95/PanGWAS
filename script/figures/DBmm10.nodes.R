#!/usr/bin/env Rscript
#Rscript DBmm10.nodes.R DBmm10.nodes.tag.txt overlap_dup_rmsk.txt non_overlap_dup_rmsk.txt
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if we have the right number of arguments
if (length(args) != 3) {
  stop("Please provide three input files:
       1. Path to nodes tag file (DBmm10.nodes.tag.txt)
       2. Path to overlap file (overlap_dup_with_rmsk.txt)
       3. Path to non-overlap file (non_overlap_dup_regions.txt)")
}

# Assign arguments to variables
nodes_file <- args[1]
overlap_file <- args[2]
non_overlap_file <- args[3]

# Read input files
data <- read.csv(nodes_file)
data_unique <- data %>%
  group_by(nodeID, tag) %>%
  slice(1) %>%  # Keep only the first occurrence of each unique combination
  ungroup()     # Ungroup to return to a normal data frame

# Now summarize by tag to get counts and percentages
data_summary <- data_unique %>%
  group_by(tag) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(percentage = count / sum(count) * 100) # 

# Create the pie chart
piechart <- ggplot(data_summary, aes(x = "", y = count, fill = tag)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "UNIQUE" = "#ae4371", 
    "DUP" = "#66c2a3",
    "UNPLACED" = "#7570b3"
  ), 
  labels = c("DUP" = "Duplicated", 
             "UNIQUE" = "Unique", 
             "UNPLACED" = "Unplaced")) +  # Custom labels for the legend
  theme_void() +
  theme(legend.position = "bottom") +
  geom_text(aes(label = sprintf("%.1f%%\n(%s)", 
                                percentage, 
                                format(count, big.mark = ","))),
            position = position_stack(vjust = 0.5)) +
  labs(title = "Distribution of pangenome nodes", fill = " ")
ggsave("/scratch/piechart_pangenome_distribution.png", piechart, width = 5, height = 4)

#####2. Total Length of Pangenome Nodes by Category
length_summary <- data_unique %>%
  group_by(tag) %>%
  summarise(total_length = sum(length))

length_summary$tag <- factor(length_summary$tag, 
                           levels = c("UNPLACED", "DUP", "UNIQUE"))  # Reversed order for stacking

# Create the plot
total_length_plot <- ggplot(length_summary, aes(x = "", y = total_length, fill = tag)) + 
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf), 
            fill = NA, color = "white", linewidth = 1) + 
  geom_bar(stat = "identity", position = "stack", width = 0.4) + 
  scale_fill_manual(
    values = c(
      "UNIQUE" = "#ae4371",
      "DUP" = "#66c2a3",
      "UNPLACED" = "#7570b3"
    ),
    labels = c(
      "UNIQUE" = "Unique",
      "DUP" = "Duplicated",
      "UNPLACED" = "Unplaced"
    ),
    breaks = c("UNIQUE", "DUP", "UNPLACED")
  ) +
  labs(x = "", y = "Length (bp)", title = " ", fill = " ") +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.1))  # Increased right-side expansion
  ) +
  theme(
    legend.position = "top",
    plot.margin = margin(5.5, 20, 5.5, 5.5),  # Increased right margin
    axis.text.y = element_text(margin = margin(r = 10)),  # Add more space for y-axis labels
    axis.text.x = element_text(margin = margin(t = 10))   # Add more space for x-axis labels
  )
ggsave("/scratch/total_length_of_pangenome_nodes.png", total_length_plot, width = 6, height = 2,dpi = 300)

#####3. Box Plot of Node Length by Type (DUP only)
data_filtered <- data_unique %>%
  filter(tag == "DUP")

# Box plot with geom_rect
boxplot_dup_length <- ggplot(data_filtered, aes(x = type, y = length, fill = type)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf), fill = NA, color = "white", linewidth = 1) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Bmm10" = "#ae4371", 
    "DB" = "#66c2a3",
    "DBmm10" = "#7570b3",
    "Dmm10" = "#1f78b4",
    "mm10" = "pink"
  )) +
  labs(x = "Type", y = "Node Length", 
       title = "Box Plot of Node Length by Type (DUP only)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("/scratch/boxplot_dup_node_length.png", boxplot_dup_length, width = 7, height = 5)

#####5. Categories of Duplicated Nodes (mantain all the occurances, not only the uniq nodes)
data_filtered_alldup <- data %>%
  filter(tag == "DUP")

data_duplicates_summary <- data_filtered_alldup %>%
  mutate(pos = as.numeric(pos)) %>%
  group_by(nodeID, chr) %>%
  summarise(
    pos_count = n_distinct(pos),
    .groups = 'keep'
  ) %>%
  group_by(nodeID) %>%
  summarise(
    different_chr = n_distinct(chr) > 1,
    has_multiple_pos = any(pos_count > 1)
  ) %>%
  mutate(
    category = case_when(
      different_chr & has_multiple_pos ~ "both",
      different_chr ~ "different chrm",
      has_multiple_pos ~ "different pos",
      TRUE ~ "same"
    )
  )

# Calculate percentage and count of each category
category_counts <- data_duplicates_summary %>%
  count(category) %>%
  mutate(Percentage = n / sum(n) * 100) 

dupcat <- ggplot(category_counts, aes(x = "", y = Percentage, fill = category))  +geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf), fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf),
            fill = "#66c2a3", alpha = 0.2) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c(
    "different pos" = "#4372ae", 
    "different chrm" = "#ae4372", 
    "both" = "#a7ae43",
    "same" = "#a7ae43"
  )) +
  labs(title = " ",
       x = "Duplicated nodes",
       y = "Percentage",
       fill = " ") +
  theme_minimal() +theme(legend.position = "top")+
  coord_flip()
ggsave("/scratch/duplicated_node_categories.png", dupcat, width = 6, height = 2)

#####6. Bar Plot of Duplicated Nodes by Length and Category
data_filtered$length <- as.numeric(data_filtered_alldup$length)
merged_data <- data_filtered_alldup %>%
  inner_join(data_duplicates_summary, by = "nodeID")

# Create length categories
merged_data$Length_Category <- cut(
  merged_data$length,
  breaks = c(-Inf, 100, 1000, 10000, Inf),
  labels = c("0-100", "100-1000", "1000-10000", ">10000")
)

# Summarize data by Length_Category and category
summary_data <- merged_data %>%
  group_by(Length_Category, category) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate total counts for each Length_Category
total_counts <- merged_data %>%
  group_by(Length_Category) %>%
  summarise(Total_Length_Category = n(), .groups = 'drop')
  
summary_data <- summary_data %>%
  left_join(total_counts, by = "Length_Category") %>%
  mutate(Percentage = Count / Total_Length_Category * 100)

options(scipen = 999)

#Bar plot showing counts of nodes by Length Category and Duplication Category
barplotcat <- ggplot(summary_data, aes(x = Length_Category, y = Count, fill = category)) +geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf), fill = NA, color = "white", linewidth = 1) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1, ymax = Inf),
            fill = "#66c2a3", alpha = 0.2)+
  geom_bar(stat = "identity", position = "dodge", width = 0.3) +  # Bar plot with dodge position
  scale_fill_manual(values = c("both" = "#a7ae43", 
                               "different chrm" = "#ae4372", 
                               "different pos" = "#4372ae",
                               "same" = "#888888")) +  # Custom fill colors for each category
  labs(title = " ",
       x = "Duplicated node length",
       y = "Count",
       fill = " ") +
  theme_minimal() +coord_flip() + scale_y_log10()+ theme(legend.position = "top")
  
##same plot but stack (didn't work, the logscale is crazy) 
barplotcatst <- ggplot(summary_data, aes(x = Length_Category, y = Count, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +  # Change to stacked bar
  scale_fill_manual(values = c("both" = "#a7ae43", 
                               "different chrm" = "#ae4372", 
                               "different pos" = "#4372ae",
                               "same" = "#888888")) +  # Custom fill colors for each category
  labs(title = " ",
       x = "Duplicated node length",
       y = "Count",
       fill = "Duplication Category") +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  scale_y_log10() +coord_flip() # Log-scale on y-axis

ggsave("/scratch/duplicated_nodes_by_length_and_category_dodge.png", barplotcat, width = 5, height = 3)
ggsave("/scratch/duplicated_nodes_by_length_and_category_stack.png", barplotcatst, width = 7, height = 5)

#####7. Distribution of Duplicated Node Lengths
dislengthnodes <- ggplot(data_filtered_alldup, aes(x = length)) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf), fill = NA, color = "white", linewidth = 1) +
  geom_histogram(binwidth = 50) +
  labs(title = "Distribution of Duplicated Node Lengths",
       x = "Duplicated node length",
       y = "Count") +
  theme_minimal() +
  scale_y_log10() + 
  xlim(0, 10000)
ggsave("/scratch/duplicated_node_length_distribution.png", dislengthnodes, width = 7, height = 5)

#####8. Repeat overlap plot (run before overlap_rmsk_dup.sh, the outputs are using here)
overlap <- read_tsv(overlap_file, col_names = c("chr", "start", "end", "repeat_class", "length", "overlap"))
non_overlap <- read_tsv(non_overlap_file, col_names = c("chr", "start", "end", "na", "length", "overlap"))
filtered_overlap <- overlap %>% filter(overlap != 0) #remove the rows where the overlap is 0
# Create length categories in your overlap and non_overlap data
filtered_overlap <- filtered_overlap %>%
  mutate(Length_Category = cut(
    length,
    breaks = c(-Inf, 100, 1000, 10000, Inf),
    labels = c("0-100", "100-1000", "1000-10000", ">10000")
  ))

non_overlap <- non_overlap %>%
  mutate(Length_Category = cut(
    length,
    breaks = c(-Inf, 100, 1000, 10000, Inf),
    labels = c("0-100", "100-1000", "1000-10000", ">10000")
  ))

# Calculate percentages for each length category
length_totals <- bind_rows(
  mutate(filtered_overlap, data_type = "overlap"),
  mutate(non_overlap, data_type = "non_overlap")
) %>%
  group_by(Length_Category) %>%
  summarise(total_count = n(), .groups = 'drop')

# For repeat families within overlap data
repeat_summary_by_length <- filtered_overlap %>%
  mutate(major_family = sapply(strsplit(as.character(repeat_class), "/"), `[`, 1)) %>%
  group_by(Length_Category, major_family) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(length_totals, by = "Length_Category") %>%
  mutate(Percentage = (count / total_count) * 100)

# For non-overlap data
non_overlap_summary <- non_overlap %>%
  group_by(Length_Category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(length_totals, by = "Length_Category") %>%
  mutate(
    Percentage = (count / total_count) * 100,
    major_family = "Non-overlap"
  )

# Combine the data
all_data <- bind_rows(
  repeat_summary_by_length,
  non_overlap_summary
)

# Create a blue-only color palette
color_palette <- c(
  "Non-overlap" = "#2BECCC",  # Keeping the original turquoise color
  "DNA" = "#000080",         # Navy blue
  "DNA?" = "#0000CD",        # Medium blue
  "LINE" = "#0066CC",        # Strong blue
  "LINE?" = "#1E90FF",       # Dodger blue
  "LTR" = "#87CEEB",         # Sky blue
  "LTR?" = "#ADD8E6",        # Light blue
  "RC" = "#4169E1",          # Royal blue
  "RC?" = "#6495ED",         # Cornflower blue
  "RNA" = "#B0E0E6",         # Powder blue
  "rRNA" = "#00688B",        # Deep sky blue
  "tRNA" = "#00BFFF",        # Deep sky blue
  "Simple_repeat" = "#87CEFA", # Light sky blue
  "SINE" = "#B0C4DE",        # Light steel blue
  "SINE?" = "#C6E2FF",       # Light sky blue
  "Satellite" = "#F0F8FF",   # Alice blue
  "snRNA" = "#4682B4",       # Steel blue
  "scRNA" = "#5F9EA0",       # Cadet blue
  "srpRNA" = "#00CED1",      # Dark turquoise
  "Low_complexity" = "#E0FFFF", # Light cyan
  "Other" = "#48D1CC",       # Medium turquoise
  "Unknown" = "#AFEEEE"      # Pale turquoise
)

# Ensure non-overlap appears first in the stacking
all_data <- all_data %>%
  mutate(major_family = factor(major_family, 
                              levels = c("Non-overlap", 
                                       sort(unique(major_family[major_family != "Non-overlap"])))))

# Create the plot with flipped coordinates
repeatov <- ggplot(all_data, aes(x = Length_Category, y = Percentage, fill = major_family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Length Category",
    y = "Percentage",
    fill = "Repeats"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.05)),
    breaks = seq(0, 100, by = 20)
  ) +
  coord_flip()
ggsave("/scratch/length_distribution_repeats.png", repeatov, width = 12, height = 8)
