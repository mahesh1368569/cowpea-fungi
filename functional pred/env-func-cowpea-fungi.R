# core
library(tidyverse)

# community ecology
library(vegan)
library(cluster)

# integration + networks
library(mixOmics)
library(WGCNA)

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

biom = import_biom("/Users/durgapurushothammaheshchinthalapudi/Documents/R-Git/cowpea-fungi/Input_files/sujan_its.biom")

metadata = import_qiime_sample_data("/Users/durgapurushothammaheshchinthalapudi/Documents/R-Git/cowpea-fungi/Input_files/metadata-its.txt")

tree = read_tree("/Users/durgapurushothammaheshchinthalapudi/Documents/R-Git/cowpea-fungi/Input_files/rooted_tree.nwk")

rep_fasta = readDNAStringSet("/Users/durgapurushothammaheshchinthalapudi/Documents/R-Git/cowpea-fungi/Input_files/fun-seq.fasta", format = "fasta")

sujan_biom = merge_phyloseq(biom,  metadata, tree)

colnames(tax_table(sujan_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(sujan_biom)

# ---- CHANGE this object name to your ITS microeco object ----
# its_mt <- <your_microeco_object>

otu <- meco_dataset$otu_table

# Ensure matrix
otu <- as.matrix(otu)

# Ensure orientation is taxa x samples
# (If your samples are rows, transpose)
if ("SampleID" %in% colnames(phys)) {
  if (all(rownames(otu) %in% phys$SampleID) && !all(colnames(otu) %in% phys$SampleID)) {
    otu <- t(otu)
  }
}

# Intersect samples
common <- intersect(colnames(otu), phys$SampleID)
length(common)

# Subset and order both identically
otu2  <- otu[, common, drop=FALSE]
phys2 <- phys[match(common, phys$SampleID), ]
stopifnot(all(colnames(otu2) == phys2$SampleID))

# prevalence filter: keep taxa present in >= 10% samples
prev <- rowSums(otu2 > 0) / ncol(otu2)
otu_filt <- otu2[prev >= 0.10, , drop=FALSE]

# optional abundance filter (tiny totals)
tot <- rowSums(otu_filt)
otu_filt <- otu_filt[tot >= 10, , drop=FALSE]

dim(otu_filt)

fungi_rel <- sweep(otu_filt, 2, colSums(otu_filt), "/")
fungi_rel[is.na(fungi_rel)] <- 0

fungi_hel <- decostand(t(fungi_rel), method="hellinger") 
# NOTE: vegan wants samples x taxa, so we transpose: t(fungi_rel) = samples x taxa

# CLR on samples x taxa
fungi_counts <- t(otu_filt)  # samples x taxa
fungi_clr <- log1p(fungi_counts)
fungi_clr <- fungi_clr - rowMeans(fungi_clr)   # simple CLR-ish with log1p (works well in practice)

# STEP -A1 Prepare predictor blocks (scaled numeric + model matrix for factors)

# scale numeric blocks
Z_enz   <- scale(phys2[, enz_vars, drop=FALSE])
Z_plant <- scale(phys2[, plant_vars, drop=FALSE])
Z_eco   <- scale(phys2[, eco_vars, drop=FALSE])
Z_soil  <- scale(phys2[, soil_vars, drop=FALSE])  # optional

# treatment design matrix (factors -> dummies)
Z_treat <- model.matrix(~ Genotype * Stage * Treatment, data = phys2)[,-1, drop=FALSE]

# STEP -A2 - Variance partioning
vp <- varpart(fungi_hel, Z_enz, Z_plant, Z_treat)
plot(vp)

#STEP- A3 - dbRDA model + permutation tests

mod <- capscale(fungi_hel ~ Z_enz + Z_plant + Z_treat, data = phys2)

anova(mod, permutations=999)             # overall
anova(mod, by="term", permutations=999)  # each block term

# Analysis B - Mantel + Partial mantel
# STEP - B1 - Build distance matrtices

D_fungi <- vegdist(fungi_hel, method="bray")  # fungi_hel is samples x taxa

D_enz   <- dist(scale(phys2[, enz_vars, drop=FALSE]))
D_plant <- dist(scale(phys2[, plant_vars, drop=FALSE]))
D_eco   <- dist(scale(phys2[, eco_vars, drop=FALSE]))

# Gower distance for mixed meta (controls)
D_meta <- daisy(phys2[, c("Genotype","Stage","Treatment")], metric="gower")

#STEP - B2 - mantel tests

mantel(D_fungi, D_enz, permutations=999)
mantel(D_fungi, D_eco, permutations=999)
mantel(D_fungi, D_plant, permutations=999)

# Partial Mantel (control drought/stage/genotype)
mantel.partial(D_fungi, D_enz, D_meta, permutations=999)
mantel.partial(D_fungi, D_eco, D_meta, permutations=999)
mantel.partial(D_fungi, D_plant, D_meta, permutations=999)
