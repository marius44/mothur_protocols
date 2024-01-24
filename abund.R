#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
#Corre con R 4.0.1
#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

source("http://bioconductor.org/biocLite.R")
biocLite('phyloseq')



##### Normalization #######

# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

##### ADONIS ###########


setwd("D:/Doctorado/BIOINFORMATICS/Amplicones_estratos_illumina/v6-v8")


# Assign variables for imported data
# ESTOS ARCHIVOS VIENEN DEL ANALISIS EN MOTHUR (SERVIDOR)
sharedfile = "stability.shared" # VIENE DE MOTHUR Y ES EL ULTIMO SHARED
#View(sharedfile)
taxfile = "stability.taxonomy" # VIENE DE MOTHUR Y ES EL ULTIMO TAXONOMY

mapfile = "metadata_isabel_actual_2022_3.csv" # ESTO ES EL METADATA QUE TIENE QUE PONER YURI DE SUS MUESTRAS


# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

# Import sample metadata
map <- read.csv(mapfile)

map <- sample_data(map)
# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

moth_merge <- merge_phyloseq(mothur_data, map)
moth_merge
colnames(tax_table(moth_merge))

colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

erie <- moth_merge %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

erie

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(erie))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
smin <- min(sample_sums(erie))
smean <- mean(sample_sums(erie))
smax <- max(sample_sums(erie))

#melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

erie_phylum <- erie %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","gray" 
)


# Plot 
ggplot(erie_phylum, aes(x = Depth.m, y = Abundance, fill = Phylum)) + 
  #  facet_grid(Depth~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    #   breaks = c("7/8", "8/4", "9/2", "10/6"),
    labels = c("0 m", "5 m", "10 m", "15 m", "20 m", "23 m", ">23 m"), 
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(), legend.position = "bottom",legend.direction = "horizontal") + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%)") +
  xlab("asdasd") +
  ggtitle("Phyla Composition of Isabel Lake \n Bacterial Communities by Sampling depth") 


#Ahora familias
erie_family <- erie %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
family_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


# Plot 
ggplot(erie_family, aes(x = Depth.m, y = Abundance, fill = Family)) + 
  #  facet_grid(Depth~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_colors) +
  scale_x_discrete(
    #   breaks = c("7/8", "8/4", "9/2", "10/6"),
    labels = c("0 m", "5 m", "10 m", "15 m", "20 m", "23 m", ">23 m"), 
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(), legend.position = "bottom",legend.direction = "horizontal") + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 2%)") +
  xlab("asdasd") +
  ggtitle("Family Composition of Isabel Lake \n Bacterial Communities by Sampling depth") 










# Scale reads to even depth 
erie_scale <- erie %>%
  scale_reads(round = "round") 

# Fix month levels in sample_data
sample_data(erie)$Depth.m <- factor(
  sample_data(erie)$Depth.m, 
  levels = c("0", "5", "10", "15", "20", "23", ">23") 
)


# Ordinate
erie_pcoa <- ordinate(
  physeq = erie, 
  method = "PCoA", 
  distance = "bray"
)

# Plot 
plot_ordination(
  physeq = erie,
  ordination = erie_pcoa,
  color = "Depth.m",
  #  shape = as.factor("pH"),
  title = "PCoA of Isabel Lake bacterial Communities"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Depth.m), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


##################
# Remove data points with missing metadata
erie_not_na <- erie %>%
  subset_samples(
    !is.na(Chloride) & 
      !is.na(Na) &
      !is.na(pH) & 
      !is.na(K) & 
      !is.na(Mg) &
      !is.na(Ca)
  )

bray_not_na <- phyloseq::distance(physeq = erie_not_na, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = erie_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ Chloride + Na + K + Mg + Salt + pH
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = erie_not_na, 
  ordination = cap_ord, 
  color = "Depth.m", 
  axes = c(1,2)
) + 
  aes()+#shape = as.factor(pH)) + 
  geom_point(aes(colour = Depth.m), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta")
  )

cap_plot
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

