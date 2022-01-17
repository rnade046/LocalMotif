#Package dependencies
packages = c("dendextend", "cluster")

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
height <- args[4]
clustMeasure <- args[5]

setwd(wd)

# Load matrix 
inputFile <- paste(projectName, "_p", pval, "_DistanceMatrix.tsv", sep="")
dm <- as.matrix(read.csv(inputFile, sep = "\t", header = F))
dm2 <- dm[,-ncol(dm)]
dm_dist <- as.dist(dm2)

# Plot dendrogram
outputFile <- paste(projectName, "_p", pval,"_", clustMeasure ,"_Dendrogram2.png", sep="")
png(outputFile, res=400, units = "in", width = 10, height = 7)
dend <- hclust(dm_dist, method = "ward.D2")%>% 
  as.dendrogram %>% 
  set("labels", NULL) %>%
  plot()%>%
  abline(h = height, lty = 2, col="black")
dev.off()

motifsFile <- paste(projectName, "_p", pval, "_MotifsInMatrix.tsv", sep="")
groupsFile <- paste(projectName, "_p", pval, "_", clustMeasure, "_group", sep="")

groups<-cutree(dend, h= height)
motifs <- as.vector(read.csv(motifsFile, header = F))
motifs$group<- groups

for (i in c(1:length(table(groups)))){
  groupX <- motifs$V1[motifs$group == i]
  fileName = paste(groupsFile, i, ".tsv", sep = "")
  write.table(groupX, fileName, sep = "\t", quote= F, col.names = F, row.names = F)
}




