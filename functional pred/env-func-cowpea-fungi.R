# core
library(tidyverse)

# community ecology
library(vegan)
library(cluster)

# integration + networks
BiocManager::install("mixOmics")

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

################################################# 111111  #################################################

#ANALYSIS A — Variance partitioning + dbRDA

#Question: how much fungal β-diversity is explained by enzymes vs plant traits vs treatment/stage/genotype?

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

################################################# 222222  #################################################

# Analysis B - Mantel + Partial mantel
#Question - Do samples that are functionally similar also have similar fungi, even after controlling treatment structure?
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

################################################# 3333  #################################################

#ANALYSIS C — DIABLO (mixOmics) multi-block signature

#Question: what fungi + enzymes + EcoPlate + plant features jointly discriminate Stage/Treatment/Group?

# Step C1 - prepare blocks and outcomes

X <- list(
  fungi  = scale(fungi_clr),
  enz    = scale(phys2[, enz_vars, drop=FALSE]),
  ecop   = scale(phys2[, eco_vars, drop=FALSE]),
  plant  = scale(phys2[, plant_vars, drop=FALSE])
)

Y <- phys2$Group   # or phys2$Stage, or interaction(phys2$Stage, phys2$Treatment)

# Step C2 - fit DIABLO

design <- matrix(0.1, ncol = length(X_use), nrow = length(X_use),
                 dimnames = list(names(X_use), names(X_use)))
diag(design) <- 0
design["fungi", c("enz","ecop","plant")] <- 0.7
design[c("enz","ecop","plant"), "fungi"] <- 0.7

diablo_fit <- block.splsda(X = X_use, Y = Y_use, ncomp = 2, design = design)


plotIndiv(diablo_fit, legend = TRUE, title = "DIABLO: samples")
circosPlot(diablo_fit, cutoff = 0.7)
network(diablo_fit, cutoff = 0.7)

# 0) make sure phys2 has rownames that are your sample IDs
# if your sample IDs are in a column like phys2$SampleID, do:
# rownames(phys2) <- phys2$SampleID

# 1) define a common set of samples present in ALL blocks + phys2
common_ids <- Reduce(intersect, c(list(rownames(phys2)), lapply(X, rownames)))

length(common_ids)              # how many samples remain
setdiff(rownames(phys2), common_ids)[1:10]   # examples dropped from phys2
lapply(X, function(m) setdiff(rownames(m), common_ids)[1:5])  # examples dropped per block

# 2) choose a single ordering (I usually follow phys2)
common_ids <- common_ids[match(common_ids, rownames(phys2))]
common_ids <- common_ids[!is.na(common_ids)]  # safety

# 3) subset/reorder phys2 and each block
phys2_use <- phys2[common_ids, , drop = FALSE]
X_use <- lapply(X, function(m) m[common_ids, , drop = FALSE])

# 4) rebuild Y from the aligned phys2
Y_use <- factor(phys2_use$Group)  # or Stage / interaction(...)
names(Y_use) <- rownames(phys2_use)

# 5) final sanity checks
stopifnot(Reduce(function(a,b) identical(rownames(a), rownames(b)), X_use))
stopifnot(all(rownames(X_use[[1]]) == names(Y_use)))

