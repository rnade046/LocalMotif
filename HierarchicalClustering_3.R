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


# wd <- "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies/corrNetTop2-400_coreTPD_p0.4_p5.41545109270352E-7/"
# projectName <- "corrNetTop2-400_coreTPD_p0.4"
# pval <- "5.41545109270352E-7"
# clustMeasure <- "ward.D2"
# height<-0.9
# condition<-"Groups_"
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

# png("MotifFamily_CoreTPD0.4_transp.png", res=400, units = "in", width = 7, height = 10, bg = "transparent")
# par(fg = 'white', col.axis = "white", col.lab = "white")
# dend <- hclust(dm_dist, method = "ward.D2") %>% 
#   as.dendrogram %>% 
#   color_branches(h=height) %>% 
#   set_labels(rep("",length(dm_dist))) %>%
#   plot(horiz = T) %>%
#   abline(v = height, lty = 2, col="white")
# dev.off()



