
### for windows
pathOS <- "C://"
### for linux
pathOS <- "/media/mirjam/OS/"

source("http://www.bioconductor.org/biocLite.R")
biocLite("gplots")

setwd(paste0(pathOS, "/Mimi/Stage_CBG/2.EXPRESSION_MODULE/"))

# Import Data and packages -------------------------------------------------------------
##package
library(isa2)
library(gplots)


# ISA functions -----------------------------------------------------------
source('expression_module/ISA_functions.R')

##data
rpkm_expression <- read.csv("Data/expression/colaus_pc_rpkm.csv", h=F)
gene_id         <- read.csv("Data/expression/colaus_pc_geneID.csv", h=F)
gene_name       <- read.csv("Data/expression/colaus_pc_geneName.csv", h=F)
pt_id           <- read.csv("Data/expression/pt.csv", h=F)
pt_id           <- paste("ID", pt_id$V1, sep = ".")
gene_info       <- cbind(gene_id, gene_name)
rownames(rpkm_expression) <- gene_id$V1
colnames(rpkm_expression) <- pt_id

# z-score normalization ---------------------------------------------------
rpkm.norm <- t(scale(t((rpkm_expression)))) # for individuals
rpkm.norm <- as.data.frame(scale(rpkm.norm)) # for rpkm values

rpkm_log <- log10(rpkm_expression+1)
rpkm_log_norm <- log10

# Quality control ---------------------------------------------------------
## density plot

plot(density(log10(rpkm_expression$ID.598)))
plot(density(log10(as.numeric(rpkm_expression[3,]))))
plot(density((rpkm.norm$ID.598)))
# 
# for(id in colnames(rpkm_expression)){
#   lines(density(log10(rpkm_expression[,id])), )
# }

plot(density(log10(rpkm_expression$ID.598)))
cl <- rainbow(555)
for(i in 1:ncol(rpkm_expression)){
  lines(density(log10(rpkm_expression[,i])), col = cl[i])
}

## normalized data
plot(density(rpkm.norm.log[,1]))
for(i in 1:ncol(rpkm.norm.log)){
  lines(density(rpkm.norm.log[,i]), col = cl[i])
}

hist(log10(as.numeric(rpkm.norm[,1])), 50)

plot(density(log10(as.numeric(rpkm_expression[1,]))))


for(gene in rownames(rpkm_expression)){
  lines(density(log10(as.numeric(rpkm_expression[gene,]))))
}


plot(density(log10(rpkm_expression$ID.598)))
plot(density((rpkm_log$ID.598)))


for(id in colnames(rpkm_expression)){
  lines(density(log10(rpkm_expression[,id])))
}

hist(log10(as.numeric(rpkm.norm[1,])), 50)

plot(density(log10(as.numeric(rpkm_expression[1,]))))

for(gene in rownames(rpkm_expression)){
  lines(density(log10(as.numeric(rpkm_expression[gene,]))))
}

## chose all gene having a rpkm value > 0 in more than 50% of the sample
qc <- NULL
for(i in 1:nrow(rpkm_expression)){
  freq <- sum(rpkm_expression[i,] == 0, na.rm = T)/ncol(rpkm_expression)
  qc <- c(qc, freq)
}

rpkm_expression_qc <- rpkm_expression[qc < 0.5,]

# log transformed data
# qc <- NULL
# for(i in 1:nrow(rpkm_log)){
#   freq <- sum(rpkm_log[i,] == 0, na.rm = T)/ncol(rpkm_log)
#   qc <- c(qc, freq)
# }
# rpkm.log.qc <- rpkm_log[qc < 0.5,]

rpkm_log_sel <- log10(rpkm_expression_qc+1)


# z-score normalization ---------------------------------------------------
rpkm.norm <- t(scale(t(rpkm_expression_qc))) # for individuals
rpkm.norm <- scale(rpkm.norm) # for rpkm values
rpkm.norm.log <- t(scale(t(rpkm.log.qc))) # for individuals
rpkm.norm.log <- scale(rpkm.norm.log) # for rpkm values

# PCA analysis ------------------------------------------------------------
pca.exp <- prcomp(rpkm.norm.log)
plot(pca.exp$rotation[,5], pca.exp$rotation[,7], pch = 19,
     main = "Principal component analysis of z-score normalized log transformed data", 
     xlab = "PC1", ylab = "PC2")

# Run simple ISA ----------------------------------------------------------------- 
# rpkm.sel.isa <- isa(as.matrix(exp_sel_var), thr.row = seq(0, 2.5, by = 0.2), thr.col = seq(0, 2.5, by = 0.2))
# rpkm.sel.isa2 <- isa(as.matrix(exp_sel_var_mean), thr.row = seq(0, 2.5, by = 0.2), thr.col = seq(0, 2.5, by = 0.2))

rpkm.isa.norm <- isa(as.matrix(rpkm.norm), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0)
rpkm.isa <- isa(as.matrix(rpkm_expression_qc), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0)
rpkm.isa.log <- isa(as.matrix(rpkm.norm.log), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0)

#### check isa using only up!!! no need to include both up and down regulated gene! 
# interested on expressed genes!not repressed!
rpkm.isa.up.up <-  isa(as.matrix(rpkm.norm), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0, 
                       direction=c("up", "up"))
rpkm.isa.up.down <-  isa(as.matrix(rpkm.norm), thr.row = seq(0.5, 3.5, by = 0.5), thr.col = 0, 
                       direction=c("updown", "up"))

save(rpkm.isa, rpkm.isa.log, rpkm.isa.norm, rpkm.isa.up, rpkm.isa.up.down, rpkm.norm.log, rpkm.norm, rpkm_expression, 
     rpkm_expression_qc, file = paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/rpkm_isa.RData"))

load(paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/rpkm_isa.RData"))

## modules information
isaModules(data = rpkm.sel.isa, type = "isa")
isaModules(data = rpkm.sel.isa2, type = "isa")
norm.mod <- isaModules(data = rpkm.isa.norm, type = "isa")
# mod  <- mod[order(mod$freq, decreasing = T),]
mod <- isaModules(data = rpkm.isa, type = "isa")
# mod2 <- mod2[order(mod2$freq, decreasing = T),]
log.mod <- isaModules(data = rpkm.isa.log, type = "isa")
up.down.mod   <- isaModules(data = rpkm.isa.up.down, type = "isa")
up.mod   <- isaModules(data = rpkm.isa.up, type = "isa")

# file for Daniel test
load(paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/modules.RData"))
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

# write.table(unique(gene_info[,2][gene_info[,1] %in% rownames(rpkm.norm)]), paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/gene_list.txt"),quote = F, row.names = F, col.names = F)

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

# lapply(ma, write, paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/test.txt"), append=T, sep = "\t", ncolumns=1000)
# write.table(ma, paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Result/test555.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

#plot
hist(mod$colGroups, col = "red")
hist(mod$rowGroups, 20, col = "blue")

## module visualization
isa2image(data = log2(rpkm.norm), type = "isa", data.isa = rpkm.isa.norm, n = 2)
isa2image(data = rpkm.norm, type = "isa", data.isa = rpkm.isa.norm, n = 2)
isa2image(data = log2(rpkm_expression_qc), type = "isa", data.isa = rpkm.isa, n = 2)
isa2image(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa, n = 2)

isa2image(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 2)

col11 <- as.matrix(isaRowNames(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 1))
col1 <- as.matrix(isaRowNames(data = log2(rpkm.norm), type = "isa", data.isa = rpkm.isa, n = 49))
col2 <- as.matrix(isaRowNames(data = log2(rpkm.norm[1:5000,]), type = "isa", data.isa = rpkm.isa, n = 2))
col2 <- as.matrix(isaRowNames(data = log2(exp_sel_var_mean), type = "isa", data.isa = rpkm.sel.isa2, n = 2))
gene_info[,2][gene_info[,1] %in% col1]
write.table(as.matrix(gene_info[,2][gene_info[,1] %in% col1]), paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/col1.txt"), quote = F, row.names = F, col.names = F)


# Effect of normalization ------------------------------------------------

identity <- matrix(NA, nrow = nrow(mod1), ncol = nrow(mod2))
for(i in 1:nrow(mod1)){
  id.u <- isaRowNames(data = rpkm.norm, type = "isa", data.isa = rpkm.isa, n = i)
  for(j in 1:nrow(mod2)){
    id.s <- isaRowNames(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa2, n = j)
    perc <- round(ifelse(length(id.u) > length(id.s), sum(id.u %in% id.s)/length(id.u), sum(id.s %in% id.u)/length(id.s)), 2)
    identity[i,j] <- perc
  }
}
rownames(identity) <- paste0("mod.u", seq(1,nrow(urine.modules)))
colnames(identity) <- paste0("mod.s", seq(1,nrow(serum.modules)))

mod57 isa
module4 isa2

isa2image(data = log2(rpkm.norm), type = "isa", data.isa = rpkm.isa, n = 57)
isa2image(data = rpkm.norm, type = "isa", data.isa = rpkm.isa, n = 57)
isa2image(data = log2(rpkm_expression_qc), type = "isa", data.isa = rpkm.isa2, n = 4)
isa2image(data = rpkm_expression_qc, type = "isa", data.isa = rpkm.isa2, n = 4)



------------------ ??--------------------

identity <- matrix(NA, nrow = length(unique(go.sel$geneSet)), ncol = length(unique(go.sel$geneSet)))
for(i in 1:length(unique(go.sel$geneSet))){
  go <- names(go.sets.hs)[unlist(lapply(go.sel$term[go.sel$geneSet == unique(go.sel$geneSet)[i]], 
                                        grep, names(go.sets.hs)))]
  for(j in 1:length(unique(go.sel$geneSet))){
    go2 <- names(go.sets.hs)[unlist(lapply(go.sel$term[go.sel$geneSet == unique(go.sel$geneSet)[j]], 
                                           grep, names(go.sets.hs)))]
    perc <- round(sum(go %in% go2)/max(length(go), length(go2)), 2)
    identity[i,j] <- perc
  }
}
rownames(identity) <- paste0("mod", unique(go.sel$geneSet+1))
colnames(identity) <- paste0("mod", unique(go.sel$geneSet+1))

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks = c(seq(0,0.2,length=100),               # for red
               seq(0.21,0.6,length=100),           # for yellow
               seq(0.61,1,length=100))             # for green
heatmap.2(identity, 
          cellnote = identity,
          notecex= 0.7,
          main = "Consistent GO terms",
          notecol="black",
          density.info="none",
          trace="none",       
          col=my_palette,
          dendrogram = "none",
          Colv="NA",
          Rowv = "NA",
          keysize = 0.5,
          key.title = "Identity",
          key.xlab = NA,
          breaks = col_breaks)



# mouse phenotype -> database


