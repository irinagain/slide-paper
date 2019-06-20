# Created March 2nd, 2017 by I.Gaynanova
# Goal: create a set of comparisons between output of JIVE and our output on BRCA data to assess the performance

#Load the original data
# Load pre-processed TCGA BRCA data
geneExp <- read.table("TCGA-BRCA/Preprocessed_GeneExp.txt", sep = "\t")
#dim(geneExp) # 645 by 348
meth <- read.table("TCGA-BRCA/Preprocessed_Meth.txt", sep = "\t")
#dim(meth) # 574 by 348
miRNA <- read.table("TCGA-BRCA/Preprocessed_miRNA.txt", sep = "\t")
#dim(miRNA) # 423 by 348
protein <- read.table("TCGA-BRCA/Preprocessed_Protein.txt", sep = "\t")
#dim(protein) # 171 by 348

# Supervised/clinical information
#################################
# Molecular subtypes
subtype <- read.table("Preprocessed_Subtype.txt", sep = "\t") # vector of length 348, from 1 to 5
# Clinical information, matched previously, happens to have subject ids.
load("Clinical_matched_BRCA.Rda")


# Concatenate the data, form pvec, and perform standardization
# Concatenate the data and form pvec
pvec <- c(nrow(geneExp), nrow(meth), nrow(miRNA), nrow(protein))
Xall <- cbind(t(geneExp), t(meth), t(miRNA), t(protein))
source("AuxillaryFunctions.R")
out_s <- standardizeX(Xall, pvec, center = T)
X <- out_s$X
svec <- out_s$svec
# Create pcum for easier indexing of datasets
pcum <- cumsum(pvec)
pcum <- c(0, pcum)
d <- 4
n <- nrow(X)
p <- ncol(X)
save(X, pvec, pcum, file = "BRCA_d4_standardized_data.Rda")

# JIVE results
load("BRCA_jive.Rda")

jive_brca$rankJ # Joint rank - 3
jive_brca$rankA # Individual ranks - 31, 28, 27, 19
rank_jive <- jive_brca$rankJ + sum(jive_brca$rankA)

# What if I refit JIVE on its structure
structure_jive <- matrix(0, rank_jive, d)
structure_jive[1:jive_brca$rankJ,] = 1
iter = jive_brca$rankJ
for (i in 1:d){
  structure_jive[(iter+ 1):(iter + jive_brca$rankA[i]),i] = 1
  iter = iter + jive_brca$rankA[i]
}
param_jive <- est_givenranks_v4(X = X, pvec = pvec, pattern = structure_jive, k_max = 1000, eps = 1e-7)
jive_refit <- param_jive

signal_jive <- matrix(0, n, p)
signal_us <- tcrossprod(param2$U, param2$V)
rank_us <- nrow(structure)
for (i in 1:d){
  signal_jive[ , (pcum[i]+1):pcum[i+1]] <- t(jive_brca$joint[[i]])/sqrt(n*p) + t(jive_brca$individual[[i]])/sqrt(n*p)
}

jive_joint <- cbind(t(jive_brca$joint[[1]]),t(jive_brca$joint[[2]]),t(jive_brca$joint[[3]]),t(jive_brca$joint[[4]]))/sqrt(n*p)
jive_individual <- cbind(t(jive_brca$individual[[1]]),t(jive_brca$individual[[2]]),t(jive_brca$individual[[3]]),t(jive_brca$individual[[4]]))/sqrt(n*p)

save(structure_jive, jive_joint, jive_individual, param_jive, file="JIVE_BRCA_d4_results.Rda")

# Our resalts on d4, using 3 folds and very samll eps
load("BRCA_d4_struct_v2.Rda")

structure = outbcv$structure_min
nrow(structure) # total rank, 50
sum(rowSums(structure) == 4) # 3, joint rank
# partial
sum(rowSums(structure) == 3) # 1, joint + partial 3 gives 4, this is strange as inconsistent with what we had for d=3 (joint 6)
# total partial
sum(rowSums(structure) == 2)  #2
# total individaul
sum((rowSums(structure) == 1)&(structure[,1]==1)) #8
sum((rowSums(structure) == 1)&(structure[,2]==1)) #5
sum((rowSums(structure) == 1)&(structure[,3]==1)) #14
sum((rowSums(structure) == 1)&(structure[,4]==1)) #17 - protein dominates everything which again is quite strange

#total individual ranks
sum(structure[,1]==1) #14
sum(structure[,2]==1) #10
sum(structure[,3]==1) #19
sum(structure[,4]==1) #20
# Again, gene expession gets the lower rank compared to miRNA/protein (meth is the lowest), and protein gets the highest rank, seems weird, in contrast to what jive does

source("SLIDEfunctions.R")
param <- est_givenranks_v4(X = X, pvec = pvec, pattern = structure, k_max = 1000, eps = 1e-7)

save(fit_seq, out_struct, outbcv, structure, param, file="SLIDE_BRCA_d4_results.Rda")

