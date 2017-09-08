
# Import package ----------------------------------------------------------
library(isa2)
library(gplots)
library(ggplot2)


# ISA functions -----------------------------------------------------------
source('C:/Mimi/Stage_CBG/2.EXPRESSION_MODULE/expression_module/ISA_functions.R')


# Grafical parameters -----------------------------------------------------
## green, yellow and red scale used for heatmap plot (used in the case of a graphical representation of the percentage of identity for example)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks = c(seq(0,0.33,length=100),             # for red
               seq(0.34,0.66,length=100),          # for yellow
               seq(0.67,1,length=100))             # for green


# Import Data -------------------------------------------------------------

## serum data
serum <- read.csv("C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/serum.nmr.focus.all.colaus1.20160830.csv", h = F)
ppm <- paste("ppm", serum[1,2:ncol(serum)], sep ="_")
tp  <- paste("id", serum[2:nrow(serum),1], sep = "_")
serum <- serum[2:nrow(serum), 2:ncol(serum)]
colnames(serum) <- ppm
rownames(serum) <- tp

## urine data
urine <- read.csv("C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/urine.nmr.focus.all.colaus1.20161205.csv", h = F)
ppm <- paste("ppm", urine[1,2:ncol(urine)], sep ="_")
tp  <- paste("id", urine[2:nrow(urine),1], sep = "_")
urine <- urine[2:nrow(urine), 2:ncol(urine)]
colnames(urine) <- ppm
rownames(urine) <- tp

rm("ppm", "tp")
## select matching id (790 individuals having serum and urine metabolomics)
serum <- serum[rownames(serum) %in% rownames(urine),]
serum <- serum[order(rownames(serum)),]
urine <- urine[rownames(urine) %in% rownames(serum),]
urine <- urine[order(rownames(urine)),]


# Data Normalization ------------------------------------------------------

## Data Normalization (z-score normalization)
# serum
serum[serum<1]=1 # to avoid negative numbers
serum_log <- log10(serum) # log-transformed data
## z-score normalization for features
serum_feature_scaled <- scale(serum_log)
## z-score normalization for individulas
serum_ind_scaled <- t(scale(t(serum_log)))
## z-score normalization first for individuals and second for the features
serum_ind_feat_scaled <- scale(serum_ind_scaled)
## z-score normaizazion first accoring to features and second for individuals
serum_feat_ind_scaled <- t(scale(t(serum_feature_scaled)))

### urine
urine[urine<1]=1 # to avoid negative numbers
urine_log <- log10(urine)# log-transformed data
## z-score normalization for features
urine_feature_scaled <- scale(urine_log)
## z-score normalization for individulas
urine_ind_scaled <- t(scale(t(urine_log)))
## z-score normalization first for individuals and second for the features
urine_ind_feat_scaled <- scale(urine_ind_scaled)
## z-score normaizazion first accoring to features and second for individuals
urine_feat_ind_scaled <- t(scale(t(urine_feature_scaled)))

### used data -> first normalized by individuals, second by features
serum <- serum_ind_feat_scaled
urine <- urine_ind_feat_scaled

#### Phenotype data:
pheno_raw <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.raw.colaus1.20161116.csv", h = F, sep = ",", stringsAsFactors = T)
pheno_transf <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.transformed.colaus1.20161116.csv", h = F, sep = ",")
pheno_names <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/trait_names.raw.colaus1.20161116.csv", h = F, sep = ",")
pheno <- pheno_raw
colnames(pheno) <- c("ID", as.character(pheno_names$V1))
rm(pheno_names, pheno_raw)

save(urine, serum, pheno, file = "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/metabo_pheno_data.Rdata")

# Look to the data!! ------------------------------------------------------
# Know the story around the data!
# Ask concrete questions!
# Always look at the data!
# Transformation if needed!

str(urine)

## Data summary and viszualization
dim(urine)
dim(serum)

table(urine<0)

plot(sort(as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2)), decreasing = T), rev(serum[1,]), type = "h")
plot(sort(as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2)), decreasing = T), rev(serum[2,]), type = "h")

summary(urine$ppm_9.28)

summary(urine)

############ visualization #################
## raw data
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/raw_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum), main = "Raw data (serum)")
dev.off()
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/raw_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine), main = "Raw data (urine)")
dev.off()

## log transformed serum data image
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/log_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_log), main = "Log transformed data (serum)")
dev.off()

## z-score normalization for features (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_scaled_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_feature_scaled), main = "z-transformation for features (serum)")
dev.off()

## z-score normalization for individulas (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_scaled_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_ind_scaled), main = "z-transformation for individuals (serum)")
dev.off()

## z-score normalization first for individuals and second for the features (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_feat_scaled_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_ind_feat_scaled), main = "z-transformation for individuals and features (serum)")
dev.off()
## z-score normalization first for features and second for the individuals (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_ind_scaled_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_feat_ind_scaled), main = "z-transformation for features and individuals (serum)")
dev.off()



### urine
## log transformed urine data image
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/log_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_log), main = "Log transformed data (urine)")
dev.off()

## z-score normalization for features (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_scaled_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_feature_scaled), main = "z-transformation for features (urine)")
dev.off()

## z-score normalization for individulas (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_scaled_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_scaled), main = "z-transformation for individuals (urine)")
dev.off()

## z-score normalization first for individuals and second for the features (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_feat_scaled_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_feat_scaled), main = "z-transformation for individuals and features (urine)")
dev.off()
## z-score normalization first for features and second for the individuals (with image)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_ind_scaled_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_feat_scaled), main = "z-transformation for features and individuals (urine)")
dev.off()


## Correlation between features
# serum
fluid = "serum"
# urine
fluid = "urine"
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_raw",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = fluid))), main = paste0("Features correlation of raw ", fluid, " data"))
dev.off()
# z-score normalized data by features
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_feat_norm",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_feature_scaled")))), main = paste0("Features correlation of features scaled ", fluid, " data"))
dev.off()
# z-score normalized data by individuals
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_ind_norm",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_ind_scaled")))), main = paste0("Features correlation of individuals scaled ", fluid, " data"))
dev.off()
# z-score normalized data first by individuals second by features
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_ind_feat_norm",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_ind_feat_scaled")))), main = paste0("Features correlation of individuals and features scaled ", fluid, " data"))
dev.off()
# z-score normalized data first by features second by indiviuals
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_feat_ind_norm",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_feat_ind_scaled")))), main = paste0("Features correlation of features and individuals scaled ", fluid, " data"))
dev.off()






# isa run ---------------------------------------------------------

# 27 07 2017 parameters
## serum
serum.isa <- isa(as.matrix(serum), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))
## urine
urine.isa <- isa(as.matrix(urine), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

## save analysis:
# save(serum.isa, urine.isa, file = "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170727.RData")
# load("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170502.RData")
# load("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170727.RData")


# urine and serum combination
uri.ser <- cbind(urine, serum)
colnames(uri.ser) <- c(paste0("u.", colnames(urine)), paste0("s.", colnames(serum)))
uri.ser.isa <- isa(as.matrix(uri.ser))



# isa unique --------------------------------------------------------------
## Increase the number of subjects per module, reduce the total number of modules (to have less, bigger modules)
serum.norm <- isa.normalize(as.matrix(serum_ind_feat_scaled))
urine.norm <- isa.normalize(as.matrix(urine_ind_feat_scaled))

## August 2017, result saved in ser_ur_isa20170827.RData (created using the server!)
#serum
serum.isa <- isa(as.matrix(serum_ind_feat_scaled), thr.row = seq(1, 2.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))
## urine
urine.isa <- isa(as.matrix(urine_ind_feat_scaled), thr.row = seq(1, 2.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

load("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_isa20170827.RData")

serum.uni <- isa.unique(serum.norm, serum.isa, method = c("cor"), 
                        ignore.div = TRUE, cor.limit = 0.7, neg.cor = TRUE, drop.zero = TRUE)

urine.uni <- isa.unique(urine.norm, urine.isa, method = c("cor"), ignore.div = TRUE, 
                        cor.limit = 0.7, neg.cor = TRUE, drop.zero = TRUE)



# Module analysis ---------------------------------------------------------

## module info
serum.modules   <- isaModules(serum.isa, type = "isa")
urine.modules   <- isaModules(urine.isa, type = "isa")
# uri.ser.modules <- isaModules(uri.ser.isa)
serum.modules.uni <- isaModules(serum.uni, type = "isa")
urine.modules.uni <- isaModules(urine.uni, type = "isa")

# serum
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(serum.modules$colGroups, col = "red", main = "Peaks by modules", xlab = "serum", cex.lab=1.6, 
     cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(serum.modules$rowGroups, col = "orange", main = "Individuals by modules", xlab = "serum", 
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()

# urine
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(urine.modules$colGroups, col = "blue", main = "Peaks by modules", xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(urine.modules$rowGroups, col = "yellow", main = "Individuals by modules", xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()

# combination of urine and serum datasets
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_uri_ser",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(uri.ser.modules$colGroups, col = "blue", main = "Peaks by modules", xlab = "Serum and urine combined", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(uri.ser.modules$rowGroups, col = "yellow", main = "Individuals by modules", xlab = "Serum and urine combined", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()

# Identity between modules ------------------------------------------------
## identity function (in ISA_functions.R file) used to calculate the percentage of identity between modules 
## either of the same biological fluid or of modules constructed using different samples (for example urine and serum)

# Overlap of individuals or ppms present in urine modules
identity.urine <- identity(data1 = urine, data2 = urine, data.isa1 = urine.isa, data.isa2 = urine.isa, 
         modules1 = urine.modules, modules2 = urine.modules, sel = 0)
identity.urine.col <- identity(data1 = urine, data2 = urine, data.isa1 = urine.isa, data.isa2 = urine.isa, 
                           modules1 = urine.modules, modules2 = urine.modules, sel = 0, Col = TRUE)
# Overlap of individuals or ppms present in serum modules
identity.serum <- identity(data1 = serum, data2 = serum, data.isa1 = serum.isa, data.isa2 = serum.isa, 
                           modules1 = serum.modules, modules2 = serum.modules, sel = 0)
identity.serum.col <- identity(data1 = serum, data2 = serum, data.isa1 = serum.isa, data.isa2 = serum.isa, 
                               modules1 = serum.modules, modules2 = serum.modules, sel = 0, Col = TRUE)

# look for similarity of modules between urine and serum 
# (selection of row and columns containings at least one value > 0.3)
identity.urine.serum <- identity(data1 = urine, data2 = serum, data.isa1 = urine.isa, data.isa2 = serum.isa, 
                                 modules1 = urine.modules, modules2 = serum.modules, sel = 0.3)
# Overlap of individuals in modules created using the combined dataset (urine-serum)
identity.uri.ser <- identity(data1 = uri.ser, data2 = uri.ser, data.isa1 = uri.ser.isa, data.isa2 = uri.ser.isa, 
                                 modules1 = uri.ser.modules, modules2 = uri.ser.modules, sel = 0)
identity.uri.ser.sel <- identity(data1 = uri.ser, data2 = uri.ser, data.isa1 = uri.ser.isa, data.isa2 = uri.ser.isa, 
                             modules1 = uri.ser.modules, modules2 = uri.ser.modules, sel = 0.8)


## Visualization
#urine
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/individuals_similarity_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine,
          main = "Urinary module having similar individuals",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/ppm_similarity_urine",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine.col,
          main = "Urinary module having similar ppm values",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

# example of similar urine modules: 4 subjects of module 149 are present in module 147 (5 subjects). 
# 21 of the 23 ppms of module 147 are present in module 149 (24 ppms)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_urine_ex",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
isa2image(urine, data.isa = urine.isa, type = "isa", n = 147)
isa2image(urine, data.isa = urine.isa, type = "isa", n = 149)
dev.off()

# serum
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/individuals_similarity_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.serum,
          main = "Serum module having similar individuals",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/ppm_similarity_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.serum.col,
          main = "Serum module having similar ppm values",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

## urine versus serum modules
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_urine_vs_serum",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine.serum, cellnote = identity.urine.serum, notecex= 0.7,
          main = "Matching individuals in serum and urine modules", margins = c(7,7), xlab = "Serum", ylab = "Urine",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5, cexRow = 1.2, cexCol = 1.2,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks)
dev.off()


head(serum.modules)
## serum -> columns
## urine -> rows of identity.urine.serum

## change row and columns names the same as serum.modules (or urine.modules)
colnames(identity.urine.serum) <- gsub(pattern = "mod", replacement = "Module", colnames(identity.urine.serum))
rownames(identity.urine.serum) <- gsub(pattern = "mod", replacement = "Module", rownames(identity.urine.serum))

## correct the identity value of two modules by the number of individuals that are present in the modules
id_es <- matrix(NA, nrow(identity.urine.serum), ncol(identity.urine.serum))
for(j in 1:ncol(identity.urine.serum)){
  n = serum.modules[colnames(identity.urine.serum)[j],][2]
  for(i in 1:nrow(identity.urine.serum)){
    m = urine.modules[rownames(identity.urine.serum)[i],][2]
    id_es[i,j] <- identity.urine.serum[i,j]/(m+n)[,1]
  }
}

colnames(id_es) <- colnames(identity.urine.serum)
rownames(id_es) <- rownames(identity.urine.serum)

heatmap.2(round(id_es/(max(id_es)), 2), cellnote = round(id_es/(max(id_es)), 2), notecex= 0.7,
          main = "Matching individuals in serum and urine modules", margins = c(7,7), xlab = "Serum", ylab = "Urine",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5, cexRow = 1.2, cexCol = 1.2,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks)
par(mfrow = c(1,2))
isa2image(serum, data.isa = serum.isa, type = "isa", n = 197)
isa2image(urine, data.isa = urine.isa, type = "isa", n = 30)


# example of similar modules
## plot modules (example of two modules with similar individuals)
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_es",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,2))
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 11, cex = 1)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 1, cex = 1)
dev.off()

png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_es2_",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,2))
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 115, cex = 1)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 257, cex = 1)
dev.off()

# Combined urine and serum dataset
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_uri_ser_joint",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.uri.ser,
          main = "Identity between modules \n created using the combined urine and serum dataset",
          notecol="black", density.info="none", trace="none", col=my_palette, dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.5, key.title = "Identity", key.xlab = NA, breaks = col_breaks,
          labRow = FALSE, labCol = FALSE)

dev.off()
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_uri_ser_joint_sel",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.uri.ser.sel, cellnote = identity.uri.ser.sel, notecex= 0.5,
          main = "Identity between modules \n created using the combined urine and serum dataset",
          notecol="black", density.info="none", trace="none", col=my_palette, dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.5, key.title = "Identity", key.xlab = NA, breaks = col_breaks)

dev.off()

par(mfrow = c(1,2))
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 124)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 24)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 43)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 181)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 83)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 64)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 224)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 193)

isaColNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 193)

-------------

# idea: identity matrix between individuals and ppms -> one minus the other -> max of absolute value
# more differences between concordance module individuals or ppms!

-------------

### correlation matrix between metabolic feature in modules created using urine and serum 
# databases having similar individuals
## module 11 in serum and module 1 in urine (example july 2017)
isaRow = serum.isa$rows[, 115] != 0
isaCol = serum.isa$columns[, 115] != 0
module.ser <- t(as.matrix(serum[isaRow, isaCol, drop=FALSE]))

isaRow = urine.isa$rows[, 257] != 0
isaCol = urine.isa$columns[, 257] != 0
module.ur <- t(as.matrix(urine[isaRow, isaCol, drop=FALSE]))

mod.ur <- module.ur[,colnames(module.ur) %in% colnames(module.ser)]
mod.ser <- module.ser[,colnames(module.ser) %in% colnames(module.ur)]

i=4
plot(as.numeric(sapply(strsplit(rownames(mod.ur), "_"), "[[", 2)), mod.ur[,i], pch = 19, col = "yellow")
points(as.numeric(sapply(strsplit(rownames(mod.ser), "_"), "[[", 2)), mod.ser[,i], pch = 19, col = "red")


ppm.u <- as.numeric(sapply(strsplit(rownames(mod.ur), "_"), "[[", 2))
ppm.s <- as.numeric(sapply(strsplit(rownames(mod.ser), "_"), "[[", 2))

i = 1
plot(mod.ser[,i], mod.ur[,i][sapply(ppm.s, function(x) which.min(abs(ppm.u-x)))])
sapply(ppm.u, function(x) which.min(abs(ppm.s-x)))

# plot(sort(mod.ur[,1]), sort(mod.ser[,1])) no sense!

intensity_correlation <- matrix(NA, nrow(mod.ser), nrow(mod.ur))
for(i in 1:nrow(mod.ser)){
  for(j in 1:nrow(mod.ur)){
    r <- cor(mod.ser[i,], mod.ur[j,], method = "pearson")
    intensity_correlation[i,j] <- round(r,2)
  }
}
rownames(intensity_correlation) <- sapply(strsplit(rownames(mod.ser), "_"), "[[", 2)
colnames(intensity_correlation) <- sapply(strsplit(rownames(mod.ur), "_"), "[[", 2)

png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/corr_matrix",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
 
heatmap.2(intensity_correlation, cellnote = intensity_correlation, notecex= 0.8, margins = c(7,7),
          main = "Intensity correlation", cexRow = 1.2, cexCol = 1.2, xlab = "ppm Urine", ylab = "ppm Serum",
          notecol="black", density.info="none", trace="none", dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.7, key.title = "Identity", key.xlab = NA)
  dev.off()
png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/corr_plot",
           Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
plot(mod.ser["ppm_2.3066",], mod.ur[ "ppm_1.0186",], pch = 19)
abline(lm(mod.ur[12,]~mod.ser[7,]), col="red", lwd = 2)
dev.off()
hist(urine[,"ppm_1.0186"],20) ## to have an idea about the intensity distribution of this spectral position
hist(serum[,"ppm_2.3066"], 20)

# PPA ---------------------------------------------------------------------
ser_ur <- list(as.matrix(t(serum)), as.matrix(t(urine)))
metabo.ppa.def <- ppa(ser_ur) # resulting in too many modules
metabo.ppa <- ppa(ser_ur, thr.row1 = seq(1.5, 2.5, by = 0.1), thr.row2 = seq(1.5, 2.5, by = 0.1), thr.col = seq(2,3, by = 0.1))

save(metabo.ppa.def, file = "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ppa20170812.RData")
load("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ppa20170812.RData")
metabo.module.info <- isaModules(data = metabo.ppa.def, type = "ppa")

hist(metabo.module.info$colGroups)
hist(metabo.module.info2$colGroups)
hist(metabo.module.info$row1Groups)

## ppa plot

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa.def,
          n = 29, name1 = "serum", name2 = "urine")

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa2,
          n = 44, name1 = "serum", name2 = "urine")


# Cross correlation -------------------------------------------------------
ccf(x = as.matrix(data.frame(ppm = es.ppm, height = es.score)), 
    y = as.numeric(as.matrix(data.frame(ppm = hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"], 
                                        height = hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"]))), type = "correlation"))


hmdb <- read.table("C://Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/hmdb.20160809.180.slop", h = T)

ppm.corr <- matrix(NA, nrow = nrow(urine.modules), ncol = length(unique(hmdb$HMDB.ID)))
for(i in 1:nrow(urine.modules)){
  ppm.u <- as.numeric(sapply(strsplit(isaColNames(data = urine, type = "isa", data.isa = urine.isa, 
                                                  n = i), "_"), "[[", 2))
  for(j in 1:length(unique(hmdb$HMDB.ID))){
    ppm.hmdb <- hmdb$ppm[hmdb$HMDB.ID  == unique(hmdb$HMDB.ID)[j]]
    ppm.tot=NULL
    for(ppm.i in hmdb$ppm[hmdb$HMDB.ID==unique(hmdb$HMDB.ID)[j]]){
      ppm  <- ppm.u[findInterval((ppm.u), c((ppm.i-0.03), (ppm.i+0.03)), rightmost.closed = T)==1]
      ppm.tot <- unique(c(ppm.tot, ppm))
    }
    perc <- length(ppm.tot)/(length(ppm.u)+length(ppm.hmdb))
    ppm.corr[i,j] <- perc
  }
}
colnames(ppm.corr) <- unique(hmdb$HMDB.ID)
rownames(ppm.corr) <- rownames(urine.modules)

ppm.max <- data.frame(module = paste0("Module", apply(ppm.corr, MARGIN = 2, FUN = which.max)), 
                      perc = round(apply(ppm.corr, MARGIN = 2, FUN = max), 2))

### "pseudospectrum" image using isa score
es <- isaColNames(urine, data.isa = urine.isa, type ="isa", n=252)
es.ppm <- as.numeric(sapply(strsplit(isaColNames(urine, data.isa = urine.isa, type = "isa", n = 252), "_"), "[[", 2))
es.score <- urine.isa$columns[colnames(urine) %in% es, 252]
colnames(urine)[138]
match(es, colnames(urine))


### "pseudospectrum" image using isa score
es <- isaColNames(serum, data.isa = serum.isa, type ="isa", n=250)
serum.isa$columns[colnames(serum) %in% es, 250]
colnames(serum)[138]
match(es, colnames(serum))

plot(serum.isa$columns[colnames(serum) %in% es, 250], type = "h")




plot(c(es.ppm, hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"]), c(es.score, -hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"]), type = "h")
plot(hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"], hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"], type = "h")



# Phenotype ---------------------------------------------------------------
pheno           <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.raw.colaus1.20161116.csv", h = F, sep = ",", stringsAsFactors = T)
pheno_transf    <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.transformed.colaus1.20161116.csv", h = F, sep = ",")
pheno_names     <- read.csv("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/trait_names.raw.colaus1.20161116.csv", h = F, sep = ",")
colnames(pheno) <- c("ID", as.character(pheno_names$V1))
colnames(pheno_transf) <- c("ID", as.character(pheno_names$V1))
pheno$ID        <- paste0("id", pheno$ID)
pheno_transf$ID <- paste0("id_", pheno_transf$ID)
rm(pheno_names)

## phenotype information
# classes extraction function:
classes <- function(data, from, to, by, variable_name, variable_name2) {
  data_df <- data.frame(data)
  data_df$variable_class <- cut(as.numeric(unlist(data[colnames(data) %in% variable_name])), seq(from, to, by))
  return(as.data.frame(table(data_df[,colnames(data_df) %in% c("variable_class", variable_name2)])))
}

# Gender distribution with age
age_classes <- classes(data = pheno, from = 30, to = 80, by = 5, variable_name = "AGE", variable_name2 = "SEX")

# Pyramidal plot
ggplot() +
  geom_col(aes(age_classes$variable_class[age_classes$SEX == 0], 
               age_classes$Freq[age_classes$SEX == 0], fill = "Sex 0", width = 0.5)) +
  geom_col(aes(age_classes$variable_class[age_classes$SEX == 1], 
               -age_classes$Freq[age_classes$SEX == 1], fill = "Sex 1", width = 0.5)) +
  scale_fill_manual(values=c("#FF0000", "#0000CC")) +
  labs(title="Age and Gender distribution", y ="Number of subjects", x = "Age classes", fill = "Gender") +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=16, face="bold"),axis.title=element_text(size=16, face="bold"))

pheno$age_class <- cut(pheno$AGE, seq(30, 80, by = 5))
boxplot(pheno$GLUC~pheno$age_class)
boxplot(pheno$SBP~pheno$age_class)
boxplot(pheno$BMI~pheno$age_class)
boxplot(pheno$TRIG~pheno$age_class)
boxplot(pheno$HDLCH~pheno$age_class)
boxplot(pheno$CHOL~pheno$age_class)
boxplot(pheno$LDLCH~pheno$age_class)
boxplot(pheno$ADTRN~pheno$age_class)

boxplot(pheno$GLUC~pheno$SEX)
boxplot(pheno$SBP~pheno$SEX)
boxplot(pheno$BMI~pheno$SEX)
boxplot(pheno$TRIG~pheno$SEX)
boxplot(pheno$HDLCH~pheno$SEX)
boxplot(pheno$CHOL~pheno$SEX)
boxplot(pheno$LDLCH~pheno$SEX)
boxplot(pheno$ADTRN~pheno$SEX)

plot(pheno$AGE, pheno$GLUC)

model <- lm(CHOL ~ SEX + AGE + BMI + PHYACT, data = pheno)
model2 <- glm(CHOL ~ SEX + AGE + BMI + PHYACT, data = pheno)
summary(lm(CHOL ~ as.factor(SEX) + AGE + BMI + as.factor(PHYACT) + as.factor(SMK), data = pheno))
summary(lm(SBP ~ SEX + AGE + BMI + PHYACT + SMK, data = pheno))
summary(lm(TRIG ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(HDLCH ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(LDLCH ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(ADTRN ~ SEX + AGE + BMI + PHYACT, data = pheno))

summary(lm(CHOL ~ BMI, data = pheno))

## Looking if some module are related to any kind of phenotype
## selecti individuals present in the metabolomics data
phen <- pheno
phen <- phen[phen$ID %in% rownames(urine),]
phen <- phen[order(phen$ID),]
phen <- unique(phen)
attach(phen)

phen <- pheno_transf[pheno_transf$ID %in% rownames(urine),]
phen <- unique(phen)

### Variable distribution
## Look at the densities distribution, at the normal qqplot and perform the shapiro-wilk test to check normality.
covariables <- c("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP")

for(trait in covariables){
  png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
             trait, Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
  plot(density(na.omit(phen[,trait])), main = paste0(trait, " Density Plot"))
  dev.off()
  png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
             trait, "qqplot", Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
  qqnorm(phen[,trait], main = paste0(trait, " Normal Q-Q Plot"))
  qqline(phen[,trait], col = 2)
  text(-2, max(phen[,trait], na.rm = T)-(max(pheno_transf[,trait], na.rm = T)/10), 
       paste0("Shapiro-Wilk normality test \n p = ", shapiro.test(phen[,trait])$p))
  dev.off()
}

## transformed data:
for(trait in covariables[c(1,2,6)]){
  png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
             trait, "tranf_", Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
  plot(density(na.omit(pheno_transf[,trait])), main = paste0("Transformed ", trait, " Density Plot"))
  dev.off()
  png(paste0("C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
             trait, "tranf_", "qqplot", Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
  qqnorm(pheno_transf[,trait], main = paste0("Transformed ", trait, " Normal Q-Q Plot"))
  qqline(pheno_transf[,trait], col = 2)
  text(-2, max(pheno_transf[,trait], na.rm = T)-(max(pheno_transf[,trait], na.rm = T)/10), 
       paste0("Shapiro-Wilk normality test \n p = ", shapiro.test(pheno_transf[,trait][1:5000])$p))
  dev.off()
}

# serum:
id.ser <- isaRowNames(serum, data.isa = serum.isa, type = "isa", n = 1)
id.ur  <- isaRowNames(urine, data.isa = urine.isa, type = "isa", n = 1)

t.test.variable <- function(phenotype.data, metabo.data, metabo.isa, variable1){
  p_val = NULL
  for(i in 1:nrow(metabo.isa$seeddata)){
    id <- isaRowNames(data = metabo.data, data.isa = metabo.isa, type = "isa", n = i)
    if(length(id)==1) {
      p <- NA
    } else {
      if(length(na.omit(phenotype.data[colnames(phenotype.data) %in% variable1][phenotype.data$ID %in% id,])) == 1) {
        p <- NA
      } else {
        p <- format(t.test(na.omit(phenotype.data[colnames(phenotype.data) %in% variable1][phenotype.data$ID %in% id,]), 
                           phenotype.data[colnames(phenotype.data) %in% variable1][!phenotype.data$ID %in% id,])$p.value, 
                    digits = 2)
      }
    }
    p_val <- c(p_val, as.numeric(p))
  }
  return(data.frame(modules = paste0("Module", 1:nrow(metabo.isa$seeddata)), p_val))  
}

age.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "AGE")
sex.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SEX")
chol.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "CHOL")
GLUC.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "GLUC")
SBP.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SBP")
BMI.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "BMI")
TRIG.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "TRIG")
HDLCH.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "HDLCH")
LDLCH.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "LDLCH")
SMK.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SMK")
PHYACT.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "PHYACT")
ADTRN.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "ADTRN")

p_val_result <- data.frame(p_val.age = age.test$p_val, p_val.sex = sex.test$p_val, p_val.chol = chol.test$p_val, 
                           p_val.GLUC = GLUC.test$p_val, p_val.SBP = SBP.test$p_val, p_val.BMI = BMI.test$p_val, 
                           p_val.TRIG = TRIG.test$p_val, p_val.HDLCH = HDLCH.test$ p_val, p_valLDLCH = LDLCH.test$p_val, 
                           p_val.SMK = SMK.test$p_val, p_val.PHYACT = PHYACT.test$p_val, p_val.ADTRN = ADTRN.test$p_val)

rownames(p_val_result) <- paste0("Module", 1:nrow(serum.modules))
p_val_result <- p_val_result[apply(p_val_result, MARGIN = 1, function(x) any(x <= 0.05/nrow(serum.modules))), ]



plot(phen$CHOL)
points(phen$CHOL[phen$ID %in% isaRowNames(serum, data.isa = serum.isa, type = 'isa', n = 93)], col = "red", pch = 19)

t.test(phen$AGE[phen$ID %in% id.ser], phen$AGE[!phen$ID %in% id.ser])
t.test(phen$AGE[phen$ID %in% id.ur], phen$AGE[!phen$ID %in% id.ur])$p.value

plot(phen$AGE, phen$GLUC, pch = 19)
points(phen$AGE[phen$ID %in% id.ser], phen$GLC[phen$ID %in% id.ser], col = "red", pch = 19)




# Phenotype - Module correaltion ------------------------------------------

# phenotype = dependent variable ("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP")
# modules = independent variables. use the score to analyse the correlation with the phenotype!
# co-variables = "AGE", "SEX", "BMI", "SMK", "PHYACT", "ADTRN" 

  
# Multiple linear regression model with phenotype ~ module's score + covariables
covariables <- c("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP")

## selection of modules having higher robustness
urine.modules.uni[order(urine.modules.uni$rob, decreasing = T),]
urine.mod.sel <- urine.modules.uni[urine.modules.uni$rob > 40,]

## calculate score for every subjects given a certain number of features
## example for module 4 of urine

ppm4 <- isaColNames(urine, data.isa = urine.uni, type = "isa", n = 4)
module4 <- urine[, ppm4]

es4 <- 1/(apply(module4, 1, sd)/min(apply(module4, 1, sd)))
summary(lm(phen[,covariables[j]] ~ es4 + SEX + AGE + BMI + PHYACT + SMK + ADTRN, data = phen))

mod4isa <- isa(module4, thr.col = 0, thr.row = 0)

linear_regression_result <- matrix(NA, ncol(urine.uni$rows), 6)
for(j in 1:length(covariables)){
  for (i in 1:ncol(urine.uni$rows)){
    p <- summary(lm(phen[,covariables[j]] ~ urine.uni$rows[,i] + SEX + AGE + BMI + PHYACT + SMK + ADTRN, data = phen))$coefficients[,4] [2]
    linear_regression_result[i,j] <- round(p,2)
  }
}

colnames(linear_regression_result) <- covariables
rownames(linear_regression_result) <- paste0("Module", 1:ncol(urine.uni$rows))

linear_regression_result <- linear_regression_result[rownames(urine.mod.sel),]


heatmap.2(linear_regression_result, cellnote = linear_regression_result, notecex= 0.7, margins = c(7,7),
          main = "Linear regression result", cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
          notecol="black", density.info="none", trace="none", dendrogram = "none", 
          col = colorRampPalette(c("red", "green"))(n = 199), breaks = c(seq(0,0.05,length=100),             # for red
                         seq(0.3,1,length=100)), Colv="NA", Rowv = "NA", 
          keysize = 0.7, key.title = "p-value", key.xlab = NA)


i = 53 # HDL and tot Chol significant
id <- isaRowNames(urine, data.isa = urine.uni, type = "isa", n = i)

plot(phen$AGE, phen$TRIG)
points(phen$AGE[phen$ID %in% id], phen$TRIG[phen$ID %in% id], col = "red", pch = 19)

i =4

## serum
serum.mod.sel <- serum.modules.uni[serum.modules.uni$rob > 30,]

linear_regression_result <- matrix(NA, ncol(serum.uni$rows), 6)
for(j in 1:length(covariables)){
  for (i in 1:ncol(serum.uni$rows)){
    p <- summary(lm(phen[,covariables[j]] ~ serum.uni$rows[,i] + SEX + AGE + BMI + PHYACT + SMK + ADTRN, data = phen))$coefficients[,4] [2]
    linear_regression_result[i,j] <- round(p,2)
  }
}


colnames(linear_regression_result) <- covariables
rownames(linear_regression_result) <- paste0("Module", 1:ncol(serum.uni$rows))

linear_regression_result <- linear_regression_result[order(serum.modules.uni$rowGroups, decreasing = T),]


heatmap.2(linear_regression_result, cellnote = linear_regression_result, notecex= 0.7, margins = c(7,7),
          main = "Linear regression result", cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
          notecol="black", density.info="none", trace="none", dendrogram = "none", 
          col = colorRampPalette(c("red", "green"))(n = 199), breaks = c(seq(0,0.05,length=100),             # for red
                                                                         seq(0.3,1,length=100)), Colv="NA", Rowv = "NA", 
          keysize = 0.7, key.title = "p-value", key.xlab = NA)

serum.modules.uni[order(serum.modules.uni$rowGroups, decreasing = T),]
