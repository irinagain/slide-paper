# Load pre-processed TCGA BRCA data
geneExp <- read.table("Preprocessed_GeneExp.txt", sep = "\t")
#dim(geneExp) # 645 by 348
meth <- read.table("Preprocessed_Meth.txt", sep = "\t")
#dim(meth) # 574 by 348
miRNA <- read.table("Preprocessed_miRNA.txt", sep = "\t")
#dim(miRNA) # 423 by 348
protein <- read.table("Preprocessed_Protein.txt", sep = "\t")
#dim(protein) # 171 by 348
subtype <- read.table("Preprocessed_Subtype.txt", sep = "\t") # vector of length 348, from 1 to 5


# # Apply JIVE
# ################################################################################
 library(r.jive)
 jive_brca <- jive(list(geneExp,meth,miRNA,protein), method = "perm", scale = T, orthIndiv = T, est = T)
# # Gives $joint list from 1 to 4, and $individual list from 1 to 4, PCA is typically done on joint
 save(jive_brca, file="BRCA_jive.Rda")

