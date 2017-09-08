#!/usr/bin/Rscript

## load library
library(ape)

## list all .phy file present
phy_file <- list.files(pattern = "*.phy$")
## number of free parameters
npar <- c(1,2,2,3,2,3,3,4,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,9,10,10,11)

best_method <- NULL
for(file in phy_file){
  phyml <- phymltest(seqfile= file, execname="./PhyML_3.0_linux64 -o lr", append = FALSE)
  AIC <- 2*(npar-phyml)
  best_method <- c(best_method, names(which.min(AIC)))
}

method_result <- data.frame(files = phy_file, best_method, tree.n = match(best_method, names(phyml)))
write.table(method_result, "best_method_result.txt", quote = F, row.names = FALSE, col.names = F)

