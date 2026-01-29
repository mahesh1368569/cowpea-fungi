library(phyloseq)
library(biomformat)
library(BiocManager)
library(file2meco)
library(MicrobiomeStat)
library(WGCNA)
library(ggtree)
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)
library(microeco)
library(ape)
library(plyr)
library(magrittr)
library(tidygraph)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(devtools)
library(psych)
library(igraph)
library(ggpubr)
library(Hmisc)
library(minpack.lm)
library(stats4)
library(EcoSimR) # load EcoSimR library
library(devEMF)
library(NST)
library(meconetcomp)
library(magrittr)
library(igraph)
library(tidyverse)
library(ape)
library(vegan)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(circlize)

## Importing files

metad <- read.table("metadata.txt",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
sampleData <- sample_data(metad)

guild_raw = read.table("otu_table.guilds.txt",header=T,sep="\t",row.names=1)

#Load reformatted FUNGuild data into R
FG <- read.table("otu_table.guilds.txt",header=T,sep="\t",row.names=1)

#Select only the column that we need
FGotus <- select(FG, -(taxonomy:Citation.Source))

FGotumat <- as(as.matrix(FGotus), "matrix")

FGOTU <- otu_table(FGotumat, taxa_are_rows = TRUE)

FGtaxmat <- select(FG, Confidence.Ranking, Taxon, Trophic.Mode, Guild, Growth.Morphology)

FGtaxmat <- as(as.matrix(FGtaxmat),"matrix")

FGTAX = tax_table(FGtaxmat)

#Creating phyloseq object
physeq = phyloseq(FGOTU,FGTAX,sampleData)

physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)

physeq.prune.nopossible = subset_taxa(physeq.prune, Confidence.Ranking="Highly Possible")

physeq.prune.nopossible = subset_taxa(physeq.prune.nopossible, Confidence.Ranking!="-")

#Create color palette
cbbPalette <- c("#009E73","#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7", "midnightblue", "lightgreen","saddlebrown", 
                "brown", "aquamarine4","lavenderblush2","snow3", "darkblue", 
                "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", 
                "black","lightskyblue", "darkgreen", "deeppink", "khaki2", 
                "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", 
                "darksalmon", "darkblue","royalblue4", "dodgerblue3", "steelblue1", 
                "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", 
                "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", 
                "darkgrey")

# Melted data for plotting
ps_data <- psmelt(physeq.prune.nopossible)

colnames(ps_data)

# Calculate total abundance by Trophic.Mode
mode_order <- ps_data %>%
  group_by(Trophic.Mode) %>%
  summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  pull(Trophic.Mode)

# Relevel the factor based on total abundance (most abundant first = bottom of stack)
ps_data$Trophic.Mode <- factor(ps_data$Trophic.Mode, levels = mode_order)

# Plot again with reordered Trophic.Mode
FUNGuildcom_ordered <- ggplot(ps_data,
                              aes(x = Treatment, y = Abundance, fill = Trophic.Mode)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  scale_fill_manual(values = cbbPalette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.8, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

# Print the plot
FUNGuildcom_ordered

ggsave("Treatment_trophic.pdf", plot = FUNGuildcom_ordered, width = 7, height = 8, dpi = 1000)
ggsave("Treatment_trophic.jpeg", plot = FUNGuildcom_ordered, width = 7, height = 8, dpi = 1000)



##################################### pathotrophs ##############################


PT <- read.table("pathotrops-guilds.txt",header=T,sep="\t",row.names=1)

#Select only the column that we need
PTotus <- select(PT, -(taxonomy:Citation.Source))

PTotumat <- as(as.matrix(PTotus), "matrix")

PTOTU <- otu_table(PTotumat, taxa_are_rows = TRUE)

PTtaxmat <- select(PT, Confidence.Ranking, Taxon, Trophic.Mode, Guild, Growth.Morphology)

PTtaxmat <- as(as.matrix(PTtaxmat),"matrix")

PTTAX = tax_table(PTtaxmat)

#Creating phyloseq object
physeq = phyloseq(PTOTU,PTTAX,sampleData)

physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)

physeq.prune.nopossible = subset_taxa(physeq.prune, Confidence.Ranking="Highly Possible")

physeq.prune.nopossible = subset_taxa(physeq.prune.nopossible, Confidence.Ranking!="-")

# Melted data for plotting
ps_data <- psmelt(physeq.prune.nopossible)

# Step 1: Summarize total abundance for each genus (Taxon)
top20_taxa <- ps_data %>%
  group_by(Taxon) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:20) %>%
  pull(Taxon)

# Step 2: Filter main dataset to include only these top 15 genera
ps_data_top15 <- ps_data %>%
  filter(Taxon %in% top20_taxa)

# Step 3: Optional - reorder factor levels for plotting
ps_data_top15$Taxon <- factor(ps_data_top15$Taxon,
                              levels = top20_taxa)

# Step 4: Plot
PT = ggplot(ps_data_top15,
            aes(x = Treatment, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  scale_fill_manual(values = cbbPalette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

PT

ggsave("patho.pdf", plot = PT, width = 8, height = 9, dpi = 1000)
ggsave("patho.jpeg", plot = PT, width = 7, height = 8, dpi = 1000)
############################################################################################################
########## circular plot visualization########################################################################


#Load reformatted FUNGuild data into R
FG <- read.table("otu_table.guilds.txt",header=T,sep="\t",row.names=1)

FGtaxmat <- select(FG, Trophic.Mode)

#Select only the column that we need
FGotus <- select(FG, -(taxonomy:Citation.Source))


# Compute relative abundance per column (sample-wise)
rel_abund <- sweep(FGotus, 2, colSums(FGotus), "/") * 100

rel_abund = cbind(rel_abund, FGtaxmat)

otu_long <- rel_abund %>%
  pivot_longer(cols = starts_with("S"), names_to = "SampleID", values_to = "Abundance") %>%
  filter(Abundance > 0)

# Load metadata with Treatment info
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t")

# Join treatment info
otu_long <- left_join(otu_long, metadata, by = c("SampleID" = "SampleID"))


mean_abundance <- otu_long %>%
  group_by(Trophic.Mode, Treatment) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

write.csv(mean_abundance, "mean-raundance-trophic.csv")

mean_tps = read.csv("mean_trophs.csv", row.names = 1, check.names = F)

data = as.matrix(mean_tps)

# Start fresh
circos.clear()

chordDiagram(data)

# Draw chord diagram with no labels first
chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 1)

# Add labels and axes
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  # Add text
  circos.text(mean(xlim), ylim[1] + 3, sector.name,# <-- distance from circle
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 1.0,  # <-- text size
              font = 2) # <-- font: 1 = plain, 2 = bold, 3 = italic
  
  # Add axis
  circos.axis(h = "top", labels.cex = 0.5,
              major.tick.percentage = 0.2,
              sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# Save directly to file before plotting
pdf('chord_plot1.pdf', width = 12, height = 12)
# Repeat the plot inside the jpeg device
circos.clear()
chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  # Add text
  circos.text(mean(xlim), ylim[1] + 2, sector.name,# <-- distance from circle
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.9,  # <-- text size
              font = 2) # <-- font: 1 = plain, 2 = bold, 3 = italic
  circos.axis(h = "top", labels.cex = 0.5,
              major.tick.percentage = 0.2,
              sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

#### Horizontal bar plot

horiplot = read.csv("mean-raundance-trophic.csv")

col_vec <- c("Control" = "red", "Cover Crop" = "#89F336")


hori = ggplot(horiplot, aes(x = MeanAbundance, 
                            y = reorder(Trophic.Mode, MeanAbundance),
                            fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = "",
       x = "Relative Abundance (%)",
       y = "Trophic Mode") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = col_vec) +  # or use scale_fill_manual() for custom colors
  theme(axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave( "horiplot.pdf", hori, height = 5, width = 8, dpi = 1000)

