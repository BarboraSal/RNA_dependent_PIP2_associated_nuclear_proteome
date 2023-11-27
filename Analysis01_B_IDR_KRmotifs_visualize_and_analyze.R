#' This code was used to perform statistical analysis and generate figures for Figure 2G-I
#' It uses the output xlsx files from the previous code as an input
#' We then summarize the results for pI, GRAVY score, and IDR length
#' Visualize the distributions and perform statistical analysis

########################################################################
# Required packages
library(openxlsx)
library(tidyverse)
library(ggridges)
library(pheatmap)

########################################################################
# define color palette
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Pastel1")

# combined result filename
analysis_label <- "Example_output"
output_file_name <- paste(Sys.Date(), analysis_label, sep = "_")

########################################################################
# Load and clean the input data
########################################################################

# setup path to working directory
setwd('c:/path/to/wd')

# get all filenames
filenames <- list.files(pattern = ".*_results\\.xlsx$")

# remove the '_results.xlsx' from filenames
all.files.names <- gsub("_results\\.xlsx$", "", filenames)

# define the substitutions as a list
substitutions <- list(
  "Cytosol_only" = "Cytosol-spec.",
  "Cytosol_all" = "Cytosolic fr.",
  "HeLa" = "Total cell prot.",
  "Human" = "Ref. human prot.",
  "Jadro_only" = "Nucleo-spec.",
  "NLIs" = "RDPA",
  "Nucleus_all" = "Nuclear fr.",
  "motif" = "M"
)

# apply the substitutions using a loop
for (original in names(substitutions)) {
  all.files.names <- gsub(original, substitutions[[original]], all.files.names)
}

# read in all the files
all.files <- lapply(filenames, function(file) {
  openxlsx::read.xlsx(file, sheet = "KR_motif_table_filt2")
})

# Set the names for the data frames in the list
names(all.files) <- all.files.names


########################################################################
# Summarize the results for pI, GRAVY score, and IDR length
########################################################################

# setup path to working directory
setwd('c:/Users/user/Documents/Barbora 2021/Projects_OLD/13_MSz_JCer_Script/_git/analysis_01/Part02_data')
list.files()

########################################################################
# combine the pI summary for visualization
pI_tab <- data.frame(matrix(ncol = 8, nrow = length(all.files)))
names(pI_tab) <- c("Dataset.motif",  "count", "pI_mean", "pI_median", "pI_SD", "pI_6","pI6_8", "pI_8")

pI_combined <- data.frame(matrix(ncol = 2, nrow = 0)) 
names(pI_combined) <- c("Dataset.motif", "pI")
 
for (i in 1:length(all.files)){
  data <- all.files[[i]]
  name <- names(all.files)[i]
  
  d <- data$pI_sequence
  
  pI_tab$Dataset.motif[i] <- name
  pI_tab$count[i] <- length(d)
  pI_tab$pI_mean[i] <- mean(d, na.rm = T)
  pI_tab$pI_median[i] <- median(d, na.rm = T)
  pI_tab$pI_SD[i] <- sd(d, na.rm = T)
  pI_tab$pI_6[i] <- sum(d < 6)
  pI_tab$pI6_8[i] <- sum(d >= 6 & d <= 8)
  pI_tab$pI_8[i] <- sum(d > 8)
  
  
  values <- data.frame(rep(name, length(d)))
  values$pI <- d
  
  names(values) <- c("Dataset.motif", "pI")
  
  pI_combined <- rbind(pI_combined, values)
  
}

########################################################################
# combine the hydrophobicity (GRAVY score) summary for visualization
gravy_tab <- data.frame(matrix(ncol = 5, nrow = length(all.files)))
names(gravy_tab) <- c("Dataset.motif",  "count", "gravy_mean", "gravy_median", "gravy_SD")

gravy_combined <- data.frame(matrix(ncol = 2, nrow = 0)) 
names(gravy_combined) <- c("Dataset.motif", "gravy")

for (i in 1:length(all.files)){
  data <- all.files[[i]]
  name <- names(all.files)[i]
  
  d <- data$hydrophobicity_sequence
  
  gravy_tab$Dataset.motif[i] <- name
  gravy_tab$count[i] <- length(d)
  gravy_tab$gravy_mean[i] <- mean(d, na.rm = T)
  gravy_tab$gravy_median[i] <- median(d, na.rm = T)
  gravy_tab$gravy_SD[i] <- sd(d, na.rm = T)
  
  
  values <- data.frame(rep(name, length(d)))
  values$gravy <- d
  
  names(values) <- c("Dataset.motif", "gravy")
  
  gravy_combined <- rbind(gravy_combined, values)
  
}

########################################################################
# combine the IDR length summary for visualization
length_tab <- data.frame(matrix(ncol = 5, nrow = length(all.files)))
names(length_tab) <- c("Dataset.motif",  "count", "length_mean", "length_median", "length_SD")

length_combined <- data.frame(matrix(ncol = 2, nrow = 0)) 
names(length_combined) <- c("Dataset.motif", "length")

for (i in 1:length(all.files)){
  data <- all.files[[i]]
  name <- names(all.files)[i]
  
  d <- data$idr_overlap_length
  
  length_tab$Dataset.motif[i] <- name
  length_tab$count[i] <- length(d)
  length_tab$length_mean[i] <- mean(d, na.rm = T)
  length_tab$length_median[i] <- median(d, na.rm = T)
  length_tab$length_SD[i] <- sd(d, na.rm = T)
  
  
  values <- data.frame(rep(name, length(d)))
  values$length <- d
  
  names(values) <- c("Dataset.motif", "length")
  
  length_combined <- rbind(length_combined, values)
  
}

########################################################################
# write to excel sheets
openxlsx::write.xlsx(x = pI_tab, file = paste(output_file_name, "pI_summary.xlsx", sep = "_"), sheetName = "pI_summary", append = FALSE)
openxlsx::write.xlsx(x = gravy_tab, file = paste(output_file_name, "GRAVY_summary.xlsx", sep = "_"), sheetName = "gravy_summary", append = FALSE)
openxlsx::write.xlsx(x = length_tab, file = paste(output_file_name, "Overlap_length_summary.xlsx", sep = "_"), sheetName = "length_summary", append = FALSE)


########################################################################
# visualize pI
########################################################################

data <- pI_combined

# reorder
names <- unique(data$Dataset.motif)

# alphabetic order
names <- names[order(names)]

# reorder based on dataset
# RDPA, nucleo-specific, cytosol-specific, total, nuclear fraction, cytosolic fraction, reference human
order <- c(7,8,9, 16, 17, 18, 13, 14, 15, 4, 5, 6, 1, 2, 3, 19, 20, 21, 10, 11, 12 )
names_order <- names[order(order)]

# create motif column for setting up different colors
data$Motif <- gsub(".*_M", "motif", data$Dataset.motif)

# plot boxplot, alphabetic
pdf(file = paste(output_file_name, "_pI_boxplot_alphabetic.pdf", sep = ""), width = 8, height = 4.5)
ggplot(data, aes(x = Dataset.motif, y = pI))+
  geom_boxplot(fill = col_pal[2])+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("pI")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# plot boxplot, reorder
pal <- col_pal[1:length(unique(data$Motif))]

pdf(file = paste(output_file_name, "_pI_boxplot_reorder.pdf", sep = ""), width = 8.5, height = 4.5)
ggplot(data, aes(x = factor(Dataset.motif, levels = names_order), y = pI, fill = Motif))+
  geom_boxplot()+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("pI")+
  scale_fill_manual(values = pal)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# pI distributions, alphabetic
pdf(file = paste(output_file_name, "_pI_ggridges_alphabetic.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(x = pI , y =  Dataset.motif, group = Dataset.motif)) + 
  geom_density_ridges(fill = col_pal[3], alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+
  xlab("pI")+
  theme(axis.title.y = element_blank())
dev.off()

# pI distributions, reorder
pdf(file = paste(output_file_name, "_pI_ggridges_reorder.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(y = factor(Dataset.motif, levels = names_order), x = pI, fill = Motif))+
  geom_density_ridges(alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+ 
  scale_fill_manual(values = pal)+
  xlab("pI")+
  theme(axis.title.y = element_blank())
dev.off()

###############################################################
# statistical analysis
test <- kruskal.test(data = data, pI ~ Dataset.motif)
test.pairwise <- pairwise.wilcox.test( x = data$pI, g = data$Dataset.motif, p.adjust.method = "BH")

write.xlsx(x = test$p.value, file = paste(output_file_name, "pI_summary.xlsx", sep = "_"), sheetName = "kruskal.pval", append = TRUE)
write.xlsx(x = test.pairwise$p.value, file = paste(output_file_name, "pI_summary.xlsx", sep = "_"), sheetName = "pairwise.wilcox.BH", append = TRUE)

###############################################################
# visualize p-values, 4 main datasets
htmp <- as.data.frame(test.pairwise$p.value)
htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_pI_htmp_4datasets.pdf", sep = "")
w = 10
h = 6

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, pI",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")

###############################################################
# visualize p-values, all datasets
htmp <- as.data.frame(test.pairwise$p.value)
# htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_pI_htmp_7datasets.pdf", sep = "")
w = 10
h = 10

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, pI",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")

########################################################################
# visualize GRAVY score
########################################################################

data <- gravy_combined
data$Motif <- gsub(".*_M", "motif", data$Dataset.motif)

# plot boxplot, alphabetic
pdf(file = paste(output_file_name, "_gravy_boxplot_alphabetic.pdf", sep = ""), width = 8, height = 4.5)
ggplot(data, aes(x = Dataset.motif, y = gravy))+
  geom_boxplot(fill = col_pal[2])+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("GRAVY")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# plot boxplot, reorder
pal <- col_pal[1:length(unique(data$Motif))]

pdf(file = paste(output_file_name, "_gravy_boxplot_reorder.pdf", sep = ""), width = 8.5, height = 4.5)
ggplot(data, aes(x = factor(Dataset.motif, levels = names_order), y = gravy, fill = Motif))+
  geom_boxplot()+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("GRAVY")+
  scale_fill_manual(values = pal)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# gravy distributions, alphabetic
pdf(file = paste(output_file_name, "_gravy_ggridges_alphabetic.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(x = gravy , y =  Dataset.motif, group = Dataset.motif)) + 
  geom_density_ridges(fill = col_pal[3], alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+
  xlab("GRAVY")+
  theme(axis.title.y = element_blank())
dev.off()

# gravy distributions, reorder
pdf(file = paste(output_file_name, "_gravy_ggridges_reorder.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(y = factor(Dataset.motif, levels = names_order), x = gravy, fill = Motif))+
  geom_density_ridges(alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+ 
  xlab("GRAVY")+
  scale_fill_manual(values = pal)+
  theme(axis.title.y = element_blank())
dev.off()

# statistical analysis
test <- kruskal.test(data = data, gravy ~ Dataset.motif)
test.pairwise <- pairwise.wilcox.test( x = data$gravy, g = data$Dataset.motif, p.adjust.method = "BH")

write.xlsx(x = test$p.value, file = paste(output_file_name, "GRAVY_summary.xlsx", sep = "_"), sheetName = "kruskal.pval", append = TRUE)
write.xlsx(x = test.pairwise$p.value, file = paste(output_file_name, "GRAVY_summary.xlsx", sep = "_"), sheetName = "pairwise.wilcox.BH", append = TRUE)

######################################################
# visualize p-values, 4 main datasets
htmp <- as.data.frame(test.pairwise$p.value)
htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_gravy_htmp_4datasets.pdf", sep = "")
w = 10
h = 6

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, GRAVY",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")

######################################################
# visualize p-values, 7 datasets
htmp <- as.data.frame(test.pairwise$p.value)
# htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_gravy_htmp_7datasets.pdf", sep = "")
w = 10
h = 10

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, GRAVY",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")

########################################################################
# visualize IDR length
########################################################################

data <- length_combined
data$Motif <- gsub(".*_M", "motif", data$Dataset.motif)

# plot boxplot, alphabetic
pdf(file = paste(output_file_name, "_length_boxplot_alphabetic.pdf", sep = ""), width = 8, height = 4.5)
ggplot(data, aes(x = Dataset.motif, y = length))+
  geom_boxplot(fill = col_pal[2])+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("IDR length")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# plot boxplot, reorder
pal <- col_pal[1:length(unique(data$Motif))]

pdf(file = paste(output_file_name, "_length_boxplot_reorder.pdf", sep = ""), width = 8.5, height = 4.5)
ggplot(data, aes(x = factor(Dataset.motif, levels = names_order), y = length, fill = Motif))+
  geom_boxplot()+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("IDR length")+
  scale_fill_manual(values = pal)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# length distributions, alphabetic
pdf(file = paste(output_file_name, "_length_ggridges_alphabetic.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(x = length , y =  Dataset.motif, group = Dataset.motif)) + 
  geom_density_ridges(fill = col_pal[3], alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+
  xlab("IDR length")+
  geom_vline(xintercept = c(250, 700))+
  theme(axis.title.y = element_blank())
dev.off()

# length distributions, reorder
pdf(file = paste(output_file_name, "_length_ggridges_reorder.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(y = factor(Dataset.motif, levels = names_order), x = length, fill = Motif))+
  geom_density_ridges(alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+ 
  xlab("IDR length")+
  scale_fill_manual(values = pal)+
  geom_vline(xintercept = c(250, 700))+
  theme(axis.title.y = element_blank())
dev.off()

######################################################
# enrichment of the long motifs in the RDPA dataset
# split based on length and count + statistical analysis
# use 250 amino acids as a threshold

data$population <- ifelse(data$length >= 250, "long", "short")

counts_total <- data %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_short <- data %>%
  filter(population == "short") %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_long <- data %>%
  filter(population == "long") %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_long <- rbind(counts_long, counts_short[1,])
counts_long <- counts_long[order(counts_long$Dataset.motif), ]
counts_long$n[1] <- 0

counts <- data.frame(Dataset.motif = counts_total$Dataset.motif, 
                     counts_total = counts_total$n, 
                     counts_short = counts_short$n, 
                     counts_long = counts_long$n)

counts$percentage_short <- (counts$counts_short/counts$counts_total)*100
counts$percentage_long <- (counts$counts_long/counts$counts_total)*100

write.xlsx(x = counts, file = paste(output_file_name, "Overlap_length_summary.xlsx", sep = "_"), sheetName = "length_count_250", append = TRUE)

######################################################
# enrichment of the long motifs in the RDPA dataset
# split based on length and count + statistical analysis
# use 700 amino acids as a threshold

data$population <- ifelse(data$length >= 700, "long", "short")

counts_total <- data %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_short <- data %>%
  filter(population == "short") %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_long <- data %>%
  filter(population == "long") %>%
  group_by(Dataset.motif) %>%
  summarise(n = n())

counts_long <- rbind(counts_long, counts_short[1:3,])
counts_long <- counts_long[order(counts_long$Dataset.motif), ]
counts_long$n[1:3] <- 0

counts <- data.frame(Dataset.motif = counts_total$Dataset.motif, 
                     counts_total = counts_total$n, 
                     counts_short = counts_short$n, 
                     counts_long = counts_long$n)

counts$percentage_short <- (counts$counts_short/counts$counts_total)*100
counts$percentage_long <- (counts$counts_long/counts$counts_total)*100

write.xlsx(x = counts, file = paste(output_file_name, "Overlap_length_summary.xlsx", sep = "_"), sheetName = "length_count_700", append = TRUE)


######################################################
# visualize log-transformed lengths
######################################################

data <- length_combined
data$Motif <- gsub(".*_M", "motif", data$Dataset.motif)

# plot boxplot, alphabetic
pdf(file = paste(output_file_name, "_log10length_boxplot_alphabetic.pdf", sep = ""), width = 8, height = 4.5)
ggplot(data, aes(x = Dataset.motif, y = log10(length)))+
  geom_boxplot(fill = col_pal[2])+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("Log10 IDR length")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# plot boxplot, reorder
pal <- col_pal[1:length(unique(data$Motif))]

pdf(file = paste(output_file_name, "_log10length_boxplot_reorder.pdf", sep = ""), width = 8.5, height = 4.5)
ggplot(data, aes(x = factor(Dataset.motif, levels = names_order), y = log10(length), fill = Motif))+
  geom_boxplot()+
  theme_classic(base_size = 14)+
  xlab("Dataset")+
  ylab("Log10 IDR length")+
  scale_fill_manual(values = pal)+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1))
dev.off()

# length distributions, alphabetic
pdf(file = paste(output_file_name, "_log10length_ggridges_alphabetic.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(x = log10(length) , y =  Dataset.motif, group = Dataset.motif)) + 
  geom_density_ridges(fill = col_pal[3], alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+
  xlab("Log10 IDR length")+
  geom_vline(xintercept = c(log10(250), log10(700)))+
  theme(axis.title.y = element_blank())
dev.off()

# length distributions, reorder
pdf(file = paste(output_file_name, "_log10length_ggridges_reorder.pdf", sep = ""), width = 6, height = 8)
ggplot(data, aes(y = factor(Dataset.motif, levels = names_order), x = log10(length), fill = Motif))+
  geom_density_ridges(alpha = 0.8, scale = 0.95)+
  theme_classic(base_size = 14)+ 
  xlab("Log10 IDR length")+
  geom_vline(xintercept = c(log10(250), log10(700)))+
  scale_fill_manual(values = pal)+
  theme(axis.title.y = element_blank())
dev.off()

# statistical analysis
str(data)

test <- kruskal.test(data = data, log10(length) ~ Dataset.motif)
test.pairwise <- pairwise.wilcox.test( x = log10(data$length), g = data$Dataset.motif, p.adjust.method = "BH")

results <- test.pairwise$p.value

write.xlsx(x = test$p.value, file = paste(output_file_name, "Overlap_length_summary.xlsx", sep = "_"), sheetName = "kruskal.pval", append = TRUE)
write.xlsx(x = test.pairwise$p.value, file = paste(output_file_name, "Overlap_length_summary.xlsx", sep = "_"), sheetName = "pairwise.wilcox.BH", append = TRUE)

########################################################################################
# visualize p-values, 4 main datasets
htmp <- as.data.frame(test.pairwise$p.value)
htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_log10length_htmp_4datasets.pdf", sep = "")
w = 10
h = 6

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, log10 length",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")

########################################################################################
# visualize p-values, 7 datasets
htmp <- as.data.frame(test.pairwise$p.value)
# htmp <- htmp[grep("spec|RDPA|Total", row.names(htmp)), grep("spec|RDPA|Total", colnames(htmp)) ]

min(htmp, na.rm = T)
max(htmp, na.rm = T)

# breakslist
breakslist <- seq(0, 0.05, by = 0.01)

# palette with breaks
cols <- c("darkred","red", "yellow", "white")
pal_rdbu_breaks <- colorRampPalette(cols)(length(breakslist))

export_name <- paste(output_file_name, "_log10length_htmp_7datasets.pdf", sep = "")
w = 10
h = 10

# heatmap
pheatmap(htmp,
         scale = "row", 
         main = "P-value matrix, log10 length",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_method = "complete", # can be the following: "ward", "single", "complete", "average", "mcquitty", "median" or       "centroid".
         color = pal_rdbu_breaks, border_color = "black",
         cellwidth = 12, cellheight = 12,
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         angle_col = c("270", "0", "45", "90", "315"), 
         display_numbers = F, number_format = "%.5f", number_color = "grey30", fontsize_number = 6,
         #annotation = annotation, 
         #annotation_colors = annotation_colors,
         #annotation_row = ann_row, 
         #annotation_col = ann_col,
         #annotation_legend = TRUE,
         #annotation_names_row = TRUE, 
         #annotation_names_col = TRUE,
         #drop_levels = TRUE, 
         show_rownames = T, show_colnames = T, 
         legend = TRUE, 
         legend_breaks = NA,  legend_labels = TRUE,
         gaps_row = 0, # gap between row colour legend
         gaps_col = 0, # gap between column colour legend
         labels_row = NULL,labels_col = NULL,
         kmeans_k = NA, breaks = breakslist, cutree_rows = NA, cutree_cols = NA,
         filename = export_name, width = w, height = h, units = "in", res = 300, 
         na_col = "#DDDDDD")
