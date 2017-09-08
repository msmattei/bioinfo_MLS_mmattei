

# Import and install ISA package ------------------------------------------------------

# source("http://bioconductor.org/biocLite.R")
# biocLite("eisa")
# biocLite("ALL")
# biocLite("biclust")
# biocLite("org.Hs.eg.db")
# install.packages("isa2")
# biocLite("pathview")
# biocLite("GO.db")

library(isa2)
# library(pathview)
# library(GO.db)


# Import Data -------------------------------------------------------------

rpkm_expression <- read.csv("Data/expression/colaus_pc_rpkm.csv", h=F)
gene_id         <- read.csv("Data/expression/colaus_pc_geneID.csv", h=F)
gene_name       <- read.csv("Data/expression/colaus_pc_geneName.csv", h=F)
pt_id           <- read.csv("Data/expression/pt.csv", h=F)
pt_id           <- paste("ID", pt_id$V1, sep = ".")
gene_info       <- cbind(gene_id, gene_name)
# colnames(rpkm_expression) <- sample_ID
# ensembl2entrez <- as.list(org.Hs.egENSEMBL2EG)
# gene_ID$entrez_id <- ensembl2entrez[gene_ID$V1]
# rownames(rpkm_expression) <- gene_ID$entrez_id
rownames(rpkm_expression) <- gene_id$V1
colnames(rpkm_expression) <- pt_id
# gene_exp <- cbind(gene_id$V1, rpkm_expression)
# write.table(gene_exp, file = "C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/expression/gene_exp20170228.txt", quote = F, sep = "\t", row.names = F)

# Quality control ---------------------------------------------------------
## chose only gene expressed at least at 10% on the overall sample --> to selective!!

# qc <- NULL
# for(i in unique(rownames(rpkm_expression))){
#   freq <- sum(rpkm_expression[i,]<1, na.rm = T)/ncol(rpkm_expression)
#   qc <- c(qc, freq)
# }
# 
# rpkm_expression_qc <- rpkm_expression[qc < 0.9,]

## chose all gene having a rpkm value < 0 in more than 50% of the sample
qc <- NULL
for(i in 1:nrow(rpkm_expression)){
  freq <- sum(rpkm_expression[i,] == 0, na.rm = T)/ncol(rpkm_expression)
  qc <- c(qc, freq)
}

rpkm_expression_qc <- rpkm_expression[qc < 0.5,]


# z-score normalization ---------------------------------------------------
rpkm.norm <- t(scale(t(rpkm_expression_qc)))
# unique(round(rowMeans(rpkm.norm), 6))
# apply(rpkm.norm, MARGIN =  1, sd)


# Calculate standard deviation, mean and variance for each gene ----------------

sd_exp        <- apply(rpkm_expression_qc, MARGIN =  1, sd) # MARGINS = 1 indicate the row (2 = columns)
mean_exp      <- apply(rpkm_expression_qc, 1, mean)
coeff_var_exp <- sd_exp/mean_exp
var_exp       <- apply(rpkm_expression_qc, 1, var)
var_mean      <- var_exp/mean_exp ## better select gene according to this value!! 
# if variance not divided by mean --> selection of most highly expressed genes!

hist(var_mean, 1000, xlim = c(1,20))
hist(log2(var_mean))

hist(coeff_var_exp, 100)
hist(log2(coeff_var_exp))

## Plot sd vs mean expression --> look for natural cutoff due to gene varaince!
smoothScatter(log2(mean_exp), log2(sd_exp), pch = 19, cex = 0.3, nrpoints = Inf)
smoothScatter(sort(log2(var_mean)), pch = 19, cex = 0.3, nrpoints = Inf)
smoothScatter(sort(log2(coeff_var_exp)), pch = 19, cex = 0.3, nrpoints = Inf)

var_top_3000      <- names(sort(coeff_var_exp, decreasing = T))[1:3000]
exp_sel_var       <- as.matrix(rpkm_expression_qc[var_top_3000,])
var_mean_top_3000 <- names(sort(var_mean, decreasing = T))[1:3000]
exp_sel_var_mean  <- as.matrix(rpkm_expression_qc[var_mean_top_3000,])

hist(log2(exp_sel_var))
smoothScatter(log2(exp_sel_var))

hist(log2(exp_sel_var_mean))
smoothScatter(log2(exp_sel_var_mean))


# PCA analysis ------------------------------------------------------------

pca.coeff.var <- prcomp(exp_sel_var)
plot(pca.coeff.var$rotation[,1], pca.coeff.var$rotation[,2], pch = 19,
     main = "Principal component analysis - gene selected by coefficient of variance", xlab = "PC1", ylab = "PC2")
pca.var.mean <- prcomp(exp_sel_var_mean)
plot(pca.var.mean$rotation[,1], pca.var.mean$rotation[,2], pch = 19,
     main = "Principal component analysis - gene selected by variance/mean", xlab = "PC1", ylab = "PC2")

pca.exp <- prcomp(rpkm.norm)
plot(pca.exp$rotation[,5], pca.exp$rotation[,7], pch = 19,
     main = "Principal component analysis of z-score transformed data", xlab = "PC1", ylab = "PC2")

# Run simple ISA ----------------------------------------------------------------- 
rpkm.sel.isa <- isa(as.matrix(exp_sel_var), thr.row = seq(0, 2.5, by = 0.2), thr.col = seq(0, 2.5, by = 0.2))
rpkm.sel.isa2 <- isa(as.matrix(exp_sel_var_mean), thr.row = seq(0, 2.5, by = 0.2), thr.col = seq(0, 2.5, by = 0.2))

rpkm.isa <- isa(as.matrix(rpkm.norm), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0)
rpkm.isa2 <- isa(as.matrix(rpkm_expression_qc), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0)

## modules information
isaModules(data = rpkm.sel.isa, type = "isa")
isaModules(data = rpkm.sel.isa2, type = "isa")
mod <- isaModules(data = rpkm.isa, type = "isa")
mod <- mod[order(mod$freq, decreasing = T),]
mod2 <- isaModules(data = rpkm.isa2, type = "isa")
mod2 <- mod2[order(mod2$freq, decreasing = T),]

# file for Daniel test
load("C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/modules.RData")
# ma <- matrix(NA, 63, 9158)
# ma[,1] <- seq(0, 62)
# ma[,2] <- rep(1, 63)
# for(i in 1:57){
#   gen <- as.matrix(isaRowNames(data = rpkm.norm, type = "isa", data.isa = rpkm.isa, n = i))
#   name <- gene_info[,2][gene_info[,1] %in% gen]
#   ma[i,3:(length(name)+2)] <- as.vector(name)
# }
# for(i in 1:6){
#   gen <- as.matrix(isaRowNames(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa2, n = i))
#   name <- gene_info[,2][gene_info[,1] %in% gen]
#   ma[(57+i),3:(length(name)+2)] <- as.vector(name)
# }

write.table(unique(gene_info[,2][gene_info[,1] %in% rownames(rpkm.norm)]), "C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/gene_list.txt",quote = F, row.names = F, col.names = F)

ma <- list()
for(i in 1:57){
  gen <- as.matrix(isaRowNames(data = rpkm.norm, type = "isa", data.isa = rpkm.isa, n = i))
  name <- gene_info[,2][gene_info[,1] %in% gen]
  ma[[i]] <- c(i-1, 1, as.vector(name))
}
for(i in 1:6){
  gen <- as.matrix(isaRowNames(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa2, n = i))
  name <- gene_info[,2][gene_info[,1] %in% gen]
  ma[[57+i]] <- c(57+i-1, 1, as.vector(name))
}

lapply(ma, write, "C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/test.txt", append=T, sep = "\t", ncolumns=1000)
# write.table(ma, "C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/test555.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#plot
hist(mod$colGroups, col = "red")
hist(mod$rowGroups, 20, col = "blue")

## module visualization
isa2image(data = log2(rpkm.norm), type = "isa", data.isa = rpkm.isa, n = 5, all = T)
isa2image(data = rpkm.norm, type = "isa", data.isa = rpkm.isa, n = 5, all = T)
isa2image(data = log2(rpkm_expression_qc), type = "isa", data.isa = rpkm.isa, n = 5, all = T)
isa2image(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa2, n = 5, all = T)

isa2image(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 2)

col11 <- as.matrix(isaRowNames(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 1))
col1 <- as.matrix(isaRowNames(data = log2(rpkm.norm), type = "isa", data.isa = rpkm.isa, n = 49))
col2 <- as.matrix(isaRowNames(data = log2(rpkm.norm[1:5000,]), type = "isa", data.isa = rpkm.isa, n = 2))
col2 <- as.matrix(isaRowNames(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 2))
gene_info[,2][gene_info[,1] %in% col1]
write.table(as.matrix(gene_info[,2][gene_info[,1] %in% col1]), "C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/col1.txt", quote = F, row.names = F, col.names = F)

