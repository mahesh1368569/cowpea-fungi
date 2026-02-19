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
library(variancePartition)


# make sure reformulas is installed (CRAN)
install.packages("reformulas")

# update Bioconductor packages (including variancePartition)
BiocManager::install("variancePartition", update = TRUE, ask = FALSE)

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

biom = import_biom("sujan_its.biom")

metadata = import_qiime_sample_data("metadata-its.txt")

tree = read_tree("rooted_tree.nwk")

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

# G+T+DS [a+b+c]:
rda.all <- rda (comm_hel ~ Genotype + Drought_Stage + Treatment, data = meta)
RsquareAdj (rda.all)

# Genotype [a+b]:
rda.Genotype <- rda (comm_hel ~ Genotype, data = meta)

RsquareAdj (rda.Genotype)

# Treatment [b+c]:
rda.Treat <- rda (comm_hel ~ Treatment, data = meta)

RsquareAdj (rda.Treat)

# Drought_Stage [b+c]:
rda.Stage <- rda (comm_hel ~ Drought_Stage, data = meta)

RsquareAdj (rda.Stage)

varp <- varpart (comm_hel, ~ Genotype, ~ Treatment, ~ Drought_Stage, data = meta)

varp

plot(varp)


anova(rda(comm_hel ~ Genotype + Condition(Treatment + Drought_Stage), data=meta))

anova(rda(comm_hel ~ Treatment + Condition(Genotype + Drought_Stage), data=meta))

anova(rda(comm_hel ~ Drought_Stage + Condition(Genotype + Treatment), data=meta))

library(vegan)
library(eulerr)

# variance partitioning
varp <- varpart(comm_hel, ~ Genotype, ~ Treatment, ~ Drought_Stage, data = meta)


ind <- varp$part$indfract  # data.frame

# helper: find the row that starts with "[a]" etc.
get_by_tag <- function(tag){
  rn <- rownames(ind)
  i <- grep(paste0("^\\[", tag, "\\]"), rn)  # tag like "a", "b", ...
  if(length(i) == 0) return(0)
  v <- ind[i[1], "Adj.R.square"]
  if(is.na(v) || !is.finite(v)) return(0)
  max(as.numeric(v), 0)  # clamp negatives to 0 for Euler
}

a <- get_by_tag("a")
b <- get_by_tag("b")
c <- get_by_tag("c")
d <- get_by_tag("d")
e <- get_by_tag("e")
f <- get_by_tag("f")
g <- get_by_tag("g")
h <- get_by_tag("h")

eul_in <- c(
  "Genotype" = a,
  "Treatment" = b,
  "Drought Stage" = c,
  "Genotype&Treatment" = d,
  "Genotype&Drought Stage" = f,
  "Treatment&Drought Stage" = e,
  "Genotype&Treatment&Drought Stage" = g
)

print(ind)                 # optional: see the table
print(round(eul_in, 6))    # should show your a,b,c as non-zero now

fit <- euler(eul_in)

pdf("varpart_euler_adjR2.pdf", width = 6.5, height = 5.5)

plot.new()

plot(
  fit,
  quantities = list(fmt = function(x) sprintf("%.3f", x)),
  fills = list(alpha = 0.5),
  edges = list(lwd = 1),
  main = ""
)
mtext(sprintf("Residual Adj.R² = %.3f", h), side = 1, line = 2)

dev.off()

##-----------------------rda-----------------------------------------##########
# fraction [a]:
rda.dose.cover <- rda (fertil.spe ~ dose + Condition (cover), data = fertil.env)
# fraction [c]:
rda.cover.dose <- rda (fertil.spe ~ cover + Condition (dose), data = fertil.env)

anova(rda.all)

anova(rda.Genotype)
anova(rda.Treat)
anova(rda.Stage)
##------mvabund model-based community test (negative bionomical deviance framework-----------##########

# comm_filt = the counts matrix you used (samples x OTUs)
# meta = metadata aligned to comm_filt rows

# drop zero-count samples
keep_samp <- rowSums(comm) > 0
comm2 <- comm[keep_samp, , drop=FALSE]
meta2 <- meta[keep_samp, , drop=FALSE]

# drop unused factor levels
meta2$Genotype <- droplevels(factor(meta2$Genotype))
meta2$Treatment  <- droplevels(factor(meta2$Treatment))
meta2$Drought_Stage    <- droplevels(factor(meta2$Drought_Stage))

# add prevalence filter (THIS is key)
prev <- colSums(comm2 > 0) / nrow(comm2)
comm2 <- comm2[, prev >= 0.10, drop=FALSE]  # try 0.15 if still unstable

# confirm no zeros
stopifnot(all(rowSums(comm2) > 0))
stopifnot(all(colSums(comm2) > 0))

library(mvabund)
Y <- mvabund(comm2)

fit_simple <- manyglm(Y ~ Genotype * Treatment * Drought_Stage,
                      data = meta2,
                      family = "negative.binomial")

anova_simple <- anova.manyglm(fit_simple,
                              p.uni = "none",
                              resamp = "montecarlo",
                              nBoot = 999)
anova_simple




# 1) Extract feature table and metadata
romm <- as(otu_table(sujan_biom), "matrix")
meta <- data.frame(phyloseq::sample_data(sujan_biom))

# make sure factors are factors
meta$Treatment     <- factor(meta$Treatment)
meta$Drought_Stage <- factor(meta$Drought_Stage)
meta$Genotype      <- factor(meta$Genotype)
meta$Replication  <- factor(meta$Replication)  

# 2) Filter very rare taxa (important for stability)
keep <- rowSums(romm > 0) >= floor(0.1 * ncol(romm))   # present in ≥10% samples
otu_f <- romm[keep, ]

# 3) CLR transform with pseudocount
otu_p <- otu_f + 0.5
rel   <- sweep(otu_p, 2, colSums(otu_p), "/")
clr   <- log(rel) - matrix(colMeans(log(rel)), nrow=nrow(rel), ncol=ncol(rel), byrow=FALSE)

dim(clr)          # features x samples
dim(meta)


head(colnames(clr))
head(rownames(meta))

length(intersect(colnames(clr), rownames(meta)))

# 4) Fit variance partition model
# Example: Treatment + Stage fixed, Plot random, Batch random
form <- ~ Treatment + Drought_Stage + Genotype

form <- ~ Treatment * Drought_Stage * Genotype + (1|Replication)

str(meta)

vp <- fitExtractVarPartModel(clr, form, meta)

vp <- sortCols(vp)

vp

plotVarPart(vp)         # genome-wide style summary (here: across taxa)
plotPercentBars(vp[1:20, ])
