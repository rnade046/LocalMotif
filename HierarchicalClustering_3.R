#Package dependencies
packages = c("dendextend", "cluster", "data.table")

package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

# Get arguments [1] = inputFile
args = commandArgs(trailingOnly=TRUE)
wd <- args[1]
projectName <- args[2]
pval <- args[3]
condition <- args[4]
height <- args[5]
clustMeasure <- args[6]


# wd <- "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies/corrNetTop2-400_TPPD_p0.3_p2.34914206101241E-8/"
# projectName <- "corrNetTop2-400_TPPD_p0.3_CoreProteins"
# pval <- "2.34914206101241E-8"
# clustMeasure <- "ward.D2"
# height<-2.7
# condition<-"Groups_CoreProteins"
setwd(wd)

# Load matrix 
inputFile <- paste(projectName, "_p", pval, "_DistanceMatrix.tsv", sep="")
dm <- as.matrix(fread(inputFile, sep = "\t", header = F))
dm2 <- dm[,-ncol(dm)]
dm_dist <- as.dist(dm2)

# Plot dendrogram
outputFile <- paste(projectName, "_p", pval,"_", clustMeasure, "_h", height ,"_Dendrogram3.png", sep="")
png(outputFile, res=400, units = "in", width = 10, height = 7)
dend <- hclust(dm_dist, method = "ward.D2") %>% 
  as.dendrogram %>% 
  set_labels(rep("",length(dm_dist))) %>%
  color_branches(h=height) %>% 
  plot()%>%
  abline(h = height, lty = 2, col="black")
dev.off()

motifsFile <- paste(projectName, "_p", pval, "_MotifsInMatrix.tsv", sep="")
groupsFile <- paste(condition, "_h", as.double(height), "/", projectName, "_p", pval, "_", clustMeasure,"_h", height, "_group", sep="")

dend <- hclust(dm_dist, method = "ward.D2")
groups<-cutree(dend, h= height)
motifs <- as.vector(read.csv(motifsFile, header = F))
motifs$group<- groups

for (i in c(1:length(table(groups)))){
  groupX <- motifs$V1[motifs$group == i]
  fileName = paste(groupsFile, i, ".tsv", sep = "")
  write.table(groupX, fileName, sep = "\t", quote= F, col.names = F, row.names = F)
}




