
# Load library ------------------------------------------------------------

library(gplots)
library(ape)

setwd("C://mmattei/UNI/Master/MLS_BIOINFORMATICS/Sequence_genome/2.IIsemester/group6HGT/")

# selection of orthologs genes ----------------------------------------
# selection of orthologs genes keeping only the one having the lowest e-value and highest percentage of alignment alignment length

output_files <- list.files("blast_output/", pattern = "*.blastout")


for (file in output_files){
  output <- read.table(paste0("blast_output/", file), h = F)
  sel.genes1 <- matrix(NA, nrow = 0, ncol = ncol(output))
  for (gene in unique(output[,1])){
    orth1 <- output[output[,1] == gene,]
    if(nrow(orth1[orth1$V11 == min(orth1$V11),]) !=1) {
      sel.genes1 <- rbind(sel.genes1, orth1[orth1$V3 == max(orth1$V3),])
    } else {
      sel.genes1 <- rbind(sel.genes1, orth1[orth1$V11 == min(orth1$V11),])
    }
  }
  sel.genes <- matrix(NA, nrow = 0, ncol = ncol(output))
  for (gene in unique(sel.genes1$V2)){
    orth <- sel.genes1[sel.genes1[,2] == gene,]
    if(nrow(orth[orth$V11 == min(orth$V11),]) !=1) {
      sel.genes <- rbind(sel.genes, orth[orth$V3 == max(orth$V3),])
    } else {
      sel.genes <- rbind(sel.genes, orth[orth$V11 == min(orth$V11),])
    }
    write.table(sel.genes, file = paste0("blast_output/gene_selection/", file), row.names = F,
                col.names = F, quote = F)
  }
}


## Calculate the ANI value for each genome comparaison
strain_list <- read.table("strain_list.txt", h = F)$V1
# strain_list <- c("Lactobacillus_bombicola", "F5_245", "F5_246", "F5_247")
output_files <- list.files("blast_output/gene_selection/", pattern = "*.blastout")

ANI <- matrix(NA, 27, 27)
colnames(ANI) <- strain_list
rownames(ANI) <- strain_list
for (file in output_files){
  output <- read.table(paste0("blast_output/gene_selection/", file), h = F)
  name <- unlist(strsplit(strsplit(file, split = "[.]")[[1]][1], "_vs_"))
  if(is.na(ANI[name[1], name[2]])){
    ANI[name[1], name[2]] <- round(mean(output$V3), 2)
  } else {
    ifelse(ANI[name[1], name[2]] > round(mean(output$V3), 2), ANI[name[1], name[2]] <- ANI[name[1], name[2]], 
           ANI[name[1], name[2]] <- round(mean(output$V3), 2))
  }
}


## heatmap of ANI value
my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 299)
col_breaks = c(seq(0,85,length=100),               
               seq(85.1,94.9,length=100),          
               seq(95,100,length=100))             
heatmap.2(ANI, cellnote = ANI, notecex= 0.7, main = "ANI value for genome comparaison", notecol="black",
          density.info="none", trace="none", dendrogram = "none",  keysize = 0.5,
          key.title = "ANI values",
          key.xlab = NA, col = my_palette, breaks = col_breaks,cexRow = 1, cexCol = 1, margins = c(6,6))



hist(ANI, 100)


## Calculate the ANI value + 2sd
strain_list <- read.table("strain_list.txt", h = F)$V1
output_files <- list.files("blast_output/gene_selection/", pattern = "*.blastout")

ANIsd <- matrix(NA, 27, 27)
colnames(ANIsd) <- strain_list
rownames(ANIsd) <- strain_list
for (file in output_files){
  output <- read.table(paste0("blast_output/gene_selection/", file), h = F)
  name <- unlist(strsplit(strsplit(file, split = "[.]")[[1]][1], "_vs_"))
  ANIsd[name[1], name[2]] <- round(mean(output$V3), 2) + round(2*(sd(output$V3)), 2)
}


# HGTransferred genes -----------------------------------------------------
output_files <- list.files("blast_output/gene_selection/", pattern = "*.blastout")

htg <- NULL
for (file in output_files){
  output <- read.table(paste0("blast_output/gene_selection/", file), h = F)
  name <- unlist(strsplit(strsplit(file, split = "[.]")[[1]][1], "_vs_"))
  hgt.gene <- output[output$V3>ANIsd[name[1], name[2]],]
  htg <- rbind(htg, hgt.gene)
}
write.table(htg, file = "results/hgt_candidate20170509.txt", row.names = F, 
            col.names = F, quote = F)



# Gene selection ----------------------------------------------------------

htg <- read.table("results/hgt_candidate20170516.txt", h = F)
orthologs <- read.table("orthologs_only", h = F)

htg_trim <- NULL
for(i in 1:length(htg$V1)){
  if(htg$V1[i] %in% unlist(orthologs) == F){
    htg_trim <- rbind(htg_trim, htg[i,])
  }
}

write.table(htg_trim, "results/hgt_candidate20170516_trim.txt", row.names = F, col.names = F, quote = F)

# Phylogeny analysis ------------------------------------------------------

## charlotte kaas file
kaas <- read.csv("../../KAAS_results.csv", h = F, sep = ";")
kaas$strain <- sapply(strsplit(x = as.character(kaas$V1), split = "|", fixed = T),"[[", 1)

out <- list()
for(id in unique(kaas$V2)){
  st <- unique(kaas$strain[kaas$V2==id])
  out[[id]] <- st
}

## Glyceraldehyde 3 phosphate DH example
setwd("C://mmattei/UNI/Master/MLS_BIOINFORMATICS/Sequence_genome/2.IIsemester/group6HGT/family/ga3pDH/")
gene <- kaas$V1[grep("K00134", kaas$V2)]
write.table(gene, file = "gad3pDH.txt", quote = F, col.names = F, row.names = F)

fam1025 <- read.table("family1025.txt", h=F)
gene %in% fam1025$V1

## run best test
test <- NULL
test <- phymltest(seqfile="family1025.phy", execname="PhyML_3.0_win32.exe -o lr", append = FALSE)
npar <- c(1,2,2,3,2,3,3,4,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,9,10,10,11)
AIC <- 2*(npar-test)
names(which.min(AIC))

## PTS system cellobiose specific
setwd("C://mmattei/UNI/Master/MLS_BIOINFORMATICS/Sequence_genome/2.IIsemester/group6HGT/family/pts/")
gene <- kaas$V1[grep("K02761", kaas$V2)]
write.table(gene, file = "/pts.txt", quote = F, col.names = F, row.names = F)

fam331 <- read.table("family331.txt", h=F)
fam573 <- read.table("family573.txt", h=F)
fam843 <- read.table("family843.txt", h=F)
fam794 <- read.table("family794.txt", h=F)
gene %in% fam331$V1
gene %in% fam573$V1 ##
gene %in% fam843$V1
gene %in% fam794$V1

## run best test
test <- NULL
test <- phymltest(seqfile="family573.phy", execname="PhyML_3.0_win32.exe -o lr", append = FALSE)
npar <- c(1,2,2,3,2,3,3,4,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,9,10,10,11)
AIC <- 2*(npar-test)
names(which.min(AIC))


## 6-phospho-beta-glucosidase example
setwd("C://mmattei/UNI/Master/MLS_BIOINFORMATICS/Sequence_genome/2.IIsemester/group6HGT/family/6pbg/")
gene <- kaas$V1[grep("K01223", kaas$V2)]
write.table(gene, file = "/pbg.txt", quote = F, col.names = F, row.names = F)

fam482 <- read.table("family482.txt", h=F)
fam1144 <- read.table("family1144.txt", h=F)
fam1036 <- read.table("family1036.txt", h=F)
gene %in% fam482$V1
gene %in% fam1144$V1
gene %in% fam1036$V1 ##

## run best test
test <- NULL
test <- phymltest(seqfile="family1036.phy", execname="PhyML_3.0_win32.exe -o lr", append = FALSE)
npar <- c(1,2,2,3,2,3,3,4,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,9,10,10,11)
AIC <- 2*(npar-test)
names(which.min(AIC))

