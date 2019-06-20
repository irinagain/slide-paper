# Load pre-processed TCGA BRCA data
geneExp <- read.table(paste(data_folder,"Preprocessed_GeneExp.txt",sep=""), sep = "\t")
#dim(geneExp) # 645 by 348
meth <- read.table(paste(data_folder,"Preprocessed_Meth.txt",sep=""), sep = "\t")
#dim(meth) # 574 by 348
miRNA <- read.table(paste(data_folder,"Preprocessed_miRNA.txt",sep=""), sep = "\t")
#dim(miRNA) # 423 by 348
protein <- read.table(paste(data_folder,"Preprocessed_Protein.txt",sep=""), sep = "\t")
#dim(protein) # 171 by 348
subtype <- read.table(paste(data_folder,"Preprocessed_Subtype.txt",sep=""), sep = "\t") # vector of length 348, from 1 to 5

# # Apply GFA
# ################################################################################

# Necessary functions for GFA
library(CCAGFA) # R package with GFA implementation
library(igraph) # for dim_select

# Concatenate the data and form pvec
pvec <- c(nrow(geneExp), nrow(meth), nrow(miRNA), nrow(protein))
Xall <- cbind(t(geneExp), t(meth), t(miRNA), t(protein))
source("../ForPaperSimulations/AuxillaryFunctions.R")

# Do standardization
out_s <- standardizeX(Xall, pvec, center = T)


# Transfer X into list form needed for GFA
# Use dim_select() for initial rank identification
######################################################
pcum <- cumsum(pvec)
data <- list()
vecr <- rep(0,4)
D = 4
for (d in 1:D){
  # Index the measurements corresponding to dth dataset
  if (d == 1){
    index <- c(1:pvec[1])
  }else{
    index <- c((pcum[d-1]+1):pcum[d])
  }
  # Allocate to corresponding list entry, columns are subjects for JIVE
  data[[d]] <- out_s$X[,index]
  vecr[d] <- dim_select(svd(data[[d]])$d)
}
norms <- out_s$norms
# vecr gives 19 6 25 12

# Apply GFA
opts = getDefaultOpts()
# Weird default choice of K
K = min(nrow(Xall), min(pvec))
out_gfa = GFAexperiment(Y = data, K = K, opts = opts, Nrep=10)
out_gfa_trim = GFAtrim(out_gfa)
# Gives K=13 total rank with everything being joint

# Set K based on total from dim_select
K = sum(vecr)
out_gfa = GFAexperiment(Y = data, K = K, opts = opts, Nrep=10)
out_gfa_trim = GFAtrim(out_gfa)
# Gives K=13 total rank with almost everything being joint

save(out_gfa, out_gfa_trim, file="BRCA_GFA.Rda")

