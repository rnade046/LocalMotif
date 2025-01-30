#Package dependencies
packages = c("dendextend", "cluster", "RColorBrewer", "data.table")

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
clustMeasure <- args[4]
lowerLimit <- args[5]
upperLimit <- args[6]
interval <- args[7]


# wd <- "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies/corrNetTop2-400_TPPD_p0.3_p2.34914206101241E-8/"
# projectName <- "corrNetTop2-400_TPPD_p0.3_CoreProteins"
# pval <- "2.34914206101241E-8"
# clustMeasure <- "ward.D2"
# 
# lowerLimit <- 2.5
# upperLimit <- 3.5
# interval <- 0.1


setwd(wd)

# Load matrix 
inputFile <- paste(projectName, "_p", pval, "_DistanceMatrix.tsv", sep="")
dm <- as.matrix(fread(inputFile, sep = "\t", header = F))
dm2 <- dm[,-ncol(dm)]
dm_dist <- as.dist(dm2)



dend <- hclust(dm_dist, method = "ward.D2")

v <- seq(as.double(lowerLimit), as.double(upperLimit), as.double(interval))
dm <- as.data.frame(v)
f <- function(x) {length(table(cutree(dend, h=x)))}
dm$groupSize <- apply(dm, 1, f)

# Plot dendrogram
outputFile <- paste(projectName, "_p", pval,"_", clustMeasure ,"_Dendrogram2_", lowerLimit, "_", upperLimit, "_", interval,".png", sep="")
png(outputFile, res=400, units = "in", width = 10, height = 7)
dend <- hclust(dm_dist, method = "ward.D2")%>% 
  as.dendrogram %>% 
  set_labels(rep("",length(dm_dist))) %>%
  plot()%>%
  abline(h = v, lty = 2, col= brewer.pal(length(v), "Set2")) %>%
  text(x = ncol(dm2) + 2, y = v, labels = dm$groupSize, cex = 0.5) %>%
  axis(side = 2, at = v, labels = F)
dev.off()
# 
# motifsFile <- paste(projectName, "_p", pval, "_MotifsInMatrix.tsv", sep="")
# groupsFile <- paste(projectName, "_p", pval, "_", clustMeasure, "_h", height, "_group", sep="")
# 
# groups<-cutree(dend, h= height)
# motifs <- as.vector(read.csv(motifsFile, header = F))
# motifs$group<- groups
# 
# for (i in c(1:length(table(groups)))){
#   groupX <- motifs$V1[motifs$group == i]
#   fileName = paste(groupsFile, i, ".tsv", sep = "")
#   write.table(groupX, fileName, sep = "\t", quote= F, col.names = F, row.names = F)
# }




