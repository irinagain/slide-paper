# Load pre-processed TCGA BRCA data
geneExp <- read.table("Preprocessed_GeneExp.txt", sep = "\t")
#dim(geneExp) # 645 by 348
meth <- read.table("Preprocessed_Meth.txt", sep = "\t")
#dim(meth) # 574 by 348
miRNA <- read.table("Preprocessed_miRNA.txt", sep = "\t")
#dim(miRNA) # 423 by 348
protein <- read.table("Preprocessed_Protein.txt", sep = "\t")
#dim(protein) # 171 by 348

# Apply our approach
#################################################
# Concatenate the data and form pvec
pvec <- c(nrow(geneExp), nrow(meth), nrow(miRNA), nrow(protein))
Xall <- cbind(t(geneExp), t(meth), t(miRNA), t(protein))
source("../AuxillaryFunctions.R")

# Do group lasso with structure first
source("../SLIDEfunctions.R")
out_s <- standardizeX(Xall, pvec, center = T)
X <- out_s$X
svec <- out_s$svec

# Form the list of candidate structures
nl = 70
fit_seq<- solve_optim1_seq_restarts(X, lambda_seq = NULL, pvec = pvec, k_max = 1000, eps = 1e-8, reduced = F, rank_total = NULL, lambda_max = max(svec), lambda_min = 0.05)

out_struct <- get_structure_v2(fit_seq, pvec)

# Select the best structure from the list
outbcv <- bcv_optim1_structure(X, pvec = pvec, structure_list = out_struct$Slist, n_fold = 3, p_fold=3, k_max = 2000, eps = 1e-8, center = T)

# Save the output
save(fit_seq, out_struct, outbcv, file="BRCA_d4_struct_v2.Rda")
