##102424 re-adjusting the figures according to reviewer requests

# Add a column as Genotype_Antibody_Sex
library(readxl)
library(stringr)
library(tidyverse)
library(ggplot2)
CD45pos_numCells <- read_excel("CD45pos_cell_types_num_cells.xlsx")
CD45pos_numCells <- CD45pos_numCells %>%
  mutate(Genotype_Antibody_Sex = stringr::str_c(Genotype, Antibody, Sex, sep="_"))
#Add columns for cluster ratio
CD45pos_numCells <- CD45pos_numCells %>%
  mutate(across(
    "B1-like":'Trem2 foamy Mac', ~.x/total, .names = "ratio_{.col}"
  ))
#Format the table for plotting
long_CD45pos_numCells <- CD45pos_numCells %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_prefix = "ratio_",
    names_to = "Cell_Type",
    values_to = "Ratio")  %>%
  group_by(Genotype_Antibody_Sex, Cell_Type) %>%
  summarize(Avg_Ratio = mean(Ratio, na.rm =TRUE), .groups = 'drop')

long_CD45pos_numCells <- long_CD45pos_numCells %>%
  separate(Genotype_Antibody_Sex, into = c("Genotype", "Antibody","Sex"), sep = "_") %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),
    Antibody = factor(Antibody, levels = c("control", "antiIl1b")),
    Genotype = factor(Genotype, levels = c("WT", "Tet2KO"))
  ) %>%
  unite("Sex_Genotype_Antibody", c("Sex", "Genotype", "Antibody"), sep="_")

View(long_CD45pos_numCells)
# Replace 'antiIl1b' with 'antiIl1β' using Unicode for beta
long_CD45pos_numCells <- long_CD45pos_numCells %>%
  mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = gsub("antiIl1b", "antiIl1\u03b2", Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))

## wrap the group names
long_CD45pos_numCells <- long_CD45pos_numCells %>%
  mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))

##re-order the groups
desired_order <- c("male_WT\ncontrol", "male_WT\nantiIl1β",
                   "female_WT\ncontrol", "female_WT\nantiIl1β",
                   "male_Tet2KO\ncontrol", "male_Tet2KO\nantiIl1β",
                   "female_Tet2KO\ncontrol", "female_Tet2KO\nantiIl1β")
long_CD45pos_numCells <- long_CD45pos_numCells %>%
  mutate(Sex_Genotype_Antibody = factor(Sex_Genotype_Antibody, levels = desired_order)) %>%
  arrange(Sex_Genotype_Antibody)




#Plotting 
CD45pos_ratio <- ggplot(long_CD45pos_numCells, aes(x=Sex_Genotype_Antibody,  y=Avg_Ratio, fill=Sex_Genotype_Antibody)) +
  geom_bar(stat="identity", position = "dodge") +
  labs(x="Sex_Genotype_Antibody", y= "Average Ratio", fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face = "bold", size=10),
        axis.text.y = element_text(size=12),  # Increases visibility of Y axis values
        axis.title.y = element_text(size=14), # Adjusts the Y axis title size
        plot.title = element_text(size=16, face="bold"),
        strip.placement = "outside", # Ensures the facet labels are positioned well
        strip.background = element_blank(), # Removes the gray background from facet labels
        strip.text.x = element_text(size=12) # Increases facet label size for visibility
        ) +
  facet_wrap(~Cell_Type, scales="free_y", ncol=3)+
  ggtitle("CD45positive_ClusterRatio")+
  scale_fill_manual(values = c("male_Tet2KO\ncontrol"="black", "male_Tet2KO\nantiIl1β"="blue", "female_Tet2KO\ncontrol"="purple", "female_Tet2KO\nantiIl1β"="magenta"))




# Read in the CD45neg cell type data
CD45neg_numCells <- read_excel("CD45neg_cell_types_num_cells.xlsx")
CD45neg_numCells <- CD45neg_numCells %>%
  mutate(Genotype_Antibody_Sex = stringr::str_c(Genotype, Antibody, Sex, sep="_"))

# Add columns for cluster ratio
CD45neg_numCells <- CD45neg_numCells %>%
  mutate(across(
    "Cardiomyocytes":'Smooth muscle cells', ~.x/total, .names = "ratio_{.col}"
  ))

# Format the table for plotting
long_CD45neg_numCells <- CD45neg_numCells %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_prefix = "ratio_",
    names_to = "Cell_Type",
    values_to = "Ratio")  %>%
  group_by(Genotype_Antibody_Sex, Cell_Type) %>%
  summarize(Avg_Ratio = mean(Ratio, na.rm =TRUE), .groups = 'drop')

long_CD45neg_numCells <- long_CD45neg_numCells %>%
  separate(Genotype_Antibody_Sex, into = c("Genotype", "Antibody","Sex"), sep = "_") %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),
    Antibody = factor(Antibody, levels = c("control", "antiIl1b")),
    Genotype = factor(Genotype, levels = c("WT", "Tet2KO"))
  ) %>%
  unite("Sex_Genotype_Antibody", c("Sex", "Genotype", "Antibody"), sep="_")

# Replace 'antiIl1b' with 'antiIl1β' using Unicode for beta
long_CD45neg_numCells <- long_CD45neg_numCells %>%
  mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = gsub("antiIl1b", "antiIl1\u03b2", Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))

# # Wrap the group names
# long_CD45neg_numCells <- long_CD45neg_numCells %>%
#   mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
#   mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))

# Re-order the groups
desired_order <- c("male_WT\ncontrol", "male_WT\nantiIl1β",
                   "female_WT\ncontrol", "female_WT\nantiIl1β",
                   "male_Tet2KO\ncontrol", "male_Tet2KO\nantiIl1β",
                   "female_Tet2KO\ncontrol", "female_Tet2KO\nantiIl1β")
long_CD45neg_numCells <- long_CD45neg_numCells %>%
  mutate(Sex_Genotype_Antibody = factor(Sex_Genotype_Antibody, levels = desired_order)) %>%
  arrange(Sex_Genotype_Antibody)

# Plotting
CD45neg_ratio <- ggplot(long_CD45neg_numCells, aes(x=Sex_Genotype_Antibody,  y=Avg_Ratio, fill=Sex_Genotype_Antibody)) +
  geom_bar(stat="identity", position = "dodge") +
  labs(x="Sex_Genotype_Antibody", y= "Average Ratio", fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face = "bold", size=10),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=16, face="bold"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12)
  ) +
  facet_wrap(~Cell_Type, scales="free_y", ncol=4) +
  ggtitle("CD45negative_ClusterRatio") +
  scale_fill_manual(values = c("male_Tet2KO\ncontrol"="black", "male_Tet2KO\nantiIl1β"="blue", "female_Tet2KO\ncontrol"="purple", "female_Tet2KO\nantiIl1β"="magenta"))

#######
##Plot the counts
##CD45 neg
# Read in the CD45neg cell type data
CD45neg_numCells <- read_excel("CD45neg_cell_types_num_cells.xlsx")
CD45neg_numCells <- CD45neg_numCells %>%
  mutate(Genotype_Antibody_Sex = stringr::str_c(Genotype, Antibody, Sex, sep="_"))

# Add a 'Total' column to the data for plotting alongside other cell types
long_CD45neg_counts <- CD45neg_numCells %>%
  pivot_longer(
    cols = c("total", "Cardiomyocytes":'Smooth muscle cells'),  # Include 'total' with other cell types
    names_to = "Cell_Type",
    values_to = "Cell_Count"
  ) %>%
  group_by(Genotype_Antibody_Sex, Cell_Type) %>%
  summarize(Avg_Cell_Count = mean(Cell_Count, na.rm = TRUE), .groups = 'drop')

# Separate columns and prepare for plotting
long_CD45neg_counts <- long_CD45neg_counts %>%
  separate(Genotype_Antibody_Sex, into = c("Genotype", "Antibody", "Sex"), sep = "_") %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),
    Antibody = factor(Antibody, levels = c("control", "antiIl1b")),
    Genotype = factor(Genotype, levels = c("WT", "Tet2KO"))
  ) %>%
  unite("Sex_Genotype_Antibody", c("Sex", "Genotype", "Antibody"), sep="_")

# Replace 'antiIl1b' with 'antiIl1β' using Unicode for beta
long_CD45neg_counts <- long_CD45neg_counts %>%
  mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = gsub("antiIl1b", "antiIl1\u03b2", Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))
write_csv(long_CD45neg_counts, "long_CD45neg_counts.csv")
# Re-order the groups
desired_order <- c("male_WT\ncontrol", "male_WT\nantiIl1β",
                   "female_WT\ncontrol", "female_WT\nantiIl1β",
                   "male_Tet2KO\ncontrol", "male_Tet2KO\nantiIl1β",
                   "female_Tet2KO\ncontrol", "female_Tet2KO\nantiIl1β")
long_CD45neg_counts <- long_CD45neg_counts %>%
  mutate(Sex_Genotype_Antibody = factor(Sex_Genotype_Antibody, levels = desired_order)) %>%
  arrange(Sex_Genotype_Antibody)

# Plotting the raw cell counts including total
CD45neg_cell_count <- ggplot(long_CD45neg_counts, aes(x=Sex_Genotype_Antibody, y=Avg_Cell_Count, fill=Sex_Genotype_Antibody)) +
  geom_bar(stat="identity", position = "dodge") +
  labs(x="Sex_Genotype_Antibody", y= "Average Cell Count", fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face = "bold", size=10),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=16, face="bold"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12)
  ) +
  facet_wrap(~Cell_Type, scales="free_y", ncol=3) +
  ggtitle("CD45negative Cell Count with Total") +
  scale_fill_manual(values = c("male_Tet2KO\ncontrol"="black", "male_Tet2KO\nantiIl1β"="blue", "female_Tet2KO\ncontrol"="purple", "female_Tet2KO\nantiIl1β"="magenta"))

##CD45pos
# Read in the CD45pos cell type data
CD45pos_numCells <- read_excel("CD45pos_cell_types_num_cells.xlsx")
CD45pos_numCells <- CD45pos_numCells %>%
  mutate(Genotype_Antibody_Sex = stringr::str_c(Genotype, Antibody, Sex, sep="_"))

# Add a 'Total' column to the data for plotting alongside other cell types
long_CD45pos_counts <- CD45pos_numCells %>%
  pivot_longer(
    cols = c("total", "B1-like":'Trem2 foamy Mac'),  # Include 'total' with other cell types
    names_to = "Cell_Type",
    values_to = "Cell_Count"
  ) %>%
  group_by(Genotype_Antibody_Sex, Cell_Type) %>%
  summarize(Avg_Cell_Count = mean(Cell_Count, na.rm = TRUE), .groups = 'drop')

# Ensure 'total' appears last by setting factor levels
cell_type_order <- c(setdiff(unique(long_CD45pos_counts$Cell_Type), "total"), "total")
long_CD45pos_counts <- long_CD45pos_counts %>%
  mutate(Cell_Type = factor(Cell_Type, levels = cell_type_order))

# Separate columns and prepare for plotting
long_CD45pos_counts <- long_CD45pos_counts %>%
  separate(Genotype_Antibody_Sex, into = c("Genotype", "Antibody", "Sex"), sep = "_") %>%
  mutate(
    Sex = factor(Sex, levels = c("female", "male")),
    Antibody = factor(Antibody, levels = c("control", "antiIl1b")),
    Genotype = factor(Genotype, levels = c("WT", "Tet2KO"))
  ) %>%
  unite("Sex_Genotype_Antibody", c("Sex", "Genotype", "Antibody"), sep="_")

# Replace 'antiIl1b' with 'antiIl1β' using Unicode for beta
long_CD45pos_counts <- long_CD45pos_counts %>%
  mutate(Sex_Genotype_Antibody = as.character(Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = gsub("antiIl1b", "antiIl1\u03b2", Sex_Genotype_Antibody)) %>%
  mutate(Sex_Genotype_Antibody = sub("_(?!.*_)", "\n", Sex_Genotype_Antibody, perl=TRUE))

# Re-order the groups
desired_order <- c("male_WT\ncontrol", "male_WT\nantiIl1β",
                   "female_WT\ncontrol", "female_WT\nantiIl1β",
                   "male_Tet2KO\ncontrol", "male_Tet2KO\nantiIl1β",
                   "female_Tet2KO\ncontrol", "female_Tet2KO\nantiIl1β")
long_CD45pos_counts <- long_CD45pos_counts %>%
  mutate(Sex_Genotype_Antibody = factor(Sex_Genotype_Antibody, levels = desired_order)) %>%
  arrange(Sex_Genotype_Antibody)

write_csv(long_CD45pos_counts, "long_CD45pos_counts.csv")
# Plotting the raw cell counts including total
CD45pos_cell_count <- ggplot(long_CD45pos_counts, aes(x=Sex_Genotype_Antibody, y=Avg_Cell_Count, fill=Sex_Genotype_Antibody)) +
  geom_bar(stat="identity", position = "dodge") +
  labs(x="Sex_Genotype_Antibody", y= "Average Cell Count", fill="Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face = "bold", size=10),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        plot.title = element_text(size=16, face="bold"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12)
  ) +
  facet_wrap(~Cell_Type, scales="free_y", ncol=4) +
  ggtitle("CD45positive Cell Count with Total") +
  scale_fill_manual(values = c("male_Tet2KO\ncontrol"="black", "male_Tet2KO\nantiIl1β"="blue", "female_Tet2KO\ncontrol"="purple", "female_Tet2KO\nantiIl1β"="magenta"))

