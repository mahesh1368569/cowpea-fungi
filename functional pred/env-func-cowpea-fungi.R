# core
library(tidyverse)

# community ecology
library(vegan)
library(cluster)

# integration + networks

library(mixOmics)
library(WGCNA)
library(vegan)
library(eulerr)
library(mvabund)



set.seed(123)

phys <- read.csv("Sujan-physiology.csv", stringsAsFactors = FALSE)

# Make sure sample IDs are unique and clean
phys$SampleID <- as.character(phys$SampleID)
stopifnot(!anyDuplicated(phys$SampleID))

# Factors (important for models)
phys$Genotype  <- factor(phys$Genotype)
phys$Stage     <- factor(phys$Stage)
phys$Treatment <- factor(phys$Treatment)

# A combined group label (useful for DIABLO / plots)
phys$Group <- interaction(phys$Genotype, phys$Stage, phys$Treatment, drop = TRUE)

# Define blocks
plant_vars <- c("gsw","E","ChlM","NFI","Plant_Height","Node_no")
enz_vars   <- c("BG","NAG")
eco_vars   <- c("AWCD","Polymers","Carbohydrates","Carboxylic_acids","Amino_acids","Amines","Phenolic_acids")
soil_vars  <- c("pH")  # optional covariate

# Build numeric matrices (scaled later)
X_plant <- phys[, plant_vars, drop=FALSE]
X_enz   <- phys[, enz_vars, drop=FALSE]
X_eco   <- phys[, eco_vars, drop=FALSE]
X_soil  <- phys[, soil_vars, drop=FALSE]

##### Importing files ########

biom = import_biom("/Users/durgasmacmini/Documents/R-Git/cowpea-fungi/Input_files/sujan_its.biom")

metadata = import_qiime_sample_data("/Users/durgasmacmini/Documents/R-Git/cowpea-fungi/Input_files/metadata-its.txt")

tree = read_tree("/Users/durgasmacmini/Documents/R-Git/cowpea-fungi/Input_files/rooted_tree.nwk")

rep_fasta = readDNAStringSet("/Users/durgasmacmini/Documents/R-Git/cowpea-fungi/Input_files/fun-seq.fasta", format = "fasta")

sujan_biom = merge_phyloseq(biom,  metadata, tree)

colnames(tax_table(sujan_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(sujan_biom)

# community matrix: samples x OTUs
comm <- t(as(otu_table(sujan_biom), "matrix"))
meta <- data.frame(phyloseq::sample_data(sujan_biom))

class(meta)

# make sure rownames match
stopifnot(all(rownames(comm) %in% rownames(meta)))
meta <- meta[rownames(comm), , drop = FALSE]

# make factors (important)
meta$Genotype    <- factor(meta$Genotype)
meta$Treatment   <- factor(meta$Treatment)
meta$Drought_Stage <- factor(meta$Drought_Stage)

comm_hel <- decostand(comm, method = "hellinger")

G <- model.matrix(~ Genotype, data = meta)[, -1, drop = FALSE]
T <- model.matrix(~ Treatment,  data = meta)[, -1, drop = FALSE]
C <- model.matrix(~ Drought_Stage,    data = meta)[, -1, drop = FALSE]

# "all interactions" set (2-way + 3-way) bundled into ONE matrix
Xfull <- model.matrix(~ Genotype * Treatment * Drought_Stage, data=meta)
Xfull <- Xfull[, colnames(Xfull) != "(Intercept)", drop=FALSE]

main_cols <- c(colnames(model.matrix(~ Genotype, data=meta)),
               colnames(model.matrix(~ Treatment, data=meta)),
               colnames(model.matrix(~ Drought_Stage, data=meta)))
main_cols <- setdiff(main_cols, "(Intercept)")

INT <- Xfull[, !colnames(Xfull) %in% main_cols, drop=FALSE]  # interactions-only

# varpart with 4 sets
vp <- varpart(comm_hel, G, T, C, INT)
vp
plot(vp)

# 1) get adjusted fractions
ind <- vp$part$indfract
fra <- ind$Adj.R.square
names(fra) <- rownames(ind)

# 2) keep only the individual fractions a–o (exclude residual p)
fra_euler <- fra[!grepl("Residuals", names(fra))]

# 3) replace negative fractions with 0 for plotting
fra_euler[fra_euler < 0] <- 0

# 4) map vegan fraction labels -> eulerr region names
# X1=Geno, X2=Treat, X3=Stage, X4=Interactions
rename_region <- function(x){
  x <- sub("^\\[[a-o]\\]\\s*=\\s*", "", x)  # remove "[a] = "
  x <- sub("\\s*\\|\\s*.*$", "", x)         # keep only left side (e.g., X1 or X1&X2)
  x <- gsub("X1", "Genotype", x)
  x <- gsub("X2", "Treatment", x)
  x <- gsub("X3", "Stage", x)
  x <- gsub("X4", "Interactions", x)
  x
}

names(fra_euler) <- vapply(names(fra_euler), rename_region, character(1))

# combine duplicates (sometimes multiple lines map to same region label)
fra_euler <- tapply(fra_euler, names(fra_euler), sum)
fra_euler <- as.numeric(fra_euler)
names(fra_euler) <- names(tapply(fra_euler, names(fra_euler), sum))  # keep names

# 5) Fit Euler
fit <- euler(fra_euler)

# 6) Plot (nice)
plot(fit,
     quantities = list(type = "percent", cex = 0.9),
     labels = list(cex = 1.1, font = 2),
     fills = list(alpha = 0.45),
     edges = list(lwd = 1.5))

# 7) Add residual text (from vp)
resid <- fra[grepl("Residuals", names(fra))]
mtext(sprintf("Residuals = %.3f", resid), side = 1, line = 1, adj = 1, cex = 0.9)

library(eulerr)

# adjusted total explained by each set (from your vp output)
# X1=G, X2=T, X3=C(Stage), X4=INT
totals <- c(
  Genotype     = vp$part$fract["[aeghklno] = X1", "Adj.R.square"],
  Treatment    = vp$part$fract["[befiklmo] = X2", "Adj.R.square"],
  Stage        = vp$part$fract["[cfgjlmno] = X3", "Adj.R.square"],
  Interactions = vp$part$fract["[dhijkmno] = X4", "Adj.R.square"]
)

totals[totals < 0] <- 0

fit_tot <- euler(totals)
plot(fit_tot,
     quantities = list(type="percent"),
     fills = list(alpha=0.45),
     edges = list(lwd=1.5))

