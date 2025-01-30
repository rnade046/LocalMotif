# Package dependencies ----------------------------------------------------
packages = c("dendextend", "cluster", "data.table", "RColorBrewer", "shiny", "miniUI")
package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = TRUE))
    }
  })

########################################################
# specify working directory file path and file names

working_directory <- "/path/to/project/motifFamilies/all/"
matrixFile <- "DistanceMatrix.tsv"
motifsFile <- "MotifsInMatrix.tsv"

########################################################

# load distance matrix ----------------------------------------------------
setwd(working_directory)
dm <- as.matrix(fread(matrixFile, sep = "\t", header = F))
dm2 <- dm[,-ncol(dm)]
dm_dist <- as.dist(dm2)
cat("Loaded matrix has : ", dim(dm2)[1], " motifs")

########################################################
# suggested clustering measures: "ward.D2", "complete", "single", "average"

clustering_measure <- "ward.D2"
########################################################

# Plot initial dendrogram -------------------------------------------------
file_name <- paste0("dendrogram1_", clustering_measure, ".png")
dend <- hclust(dm_dist, method = clustering_measure)%>% 
  as.dendrogram %>% 
  set_labels(rep("",length(dm_dist))) %>%
  plot()

# plot is output as png by default - other possible types = pdf, svg
if (file_name != "") {
  dev.copy(png, filename = file_name, res=400, units = "in", width = 10, height = 7)
  dev.off()
  cat("Plot saved as", file_name, "\n")
}

###################################################
## users should specify thresholds to evaluate based on previous plot ###

lowerLimit <- 0.5
upperLimit <- 1.5
interval <- 0.25

##################################################

dend <- hclust(dm_dist, method = clustering_measure)

v <- seq(as.double(lowerLimit), as.double(upperLimit), as.double(interval))
dm <- as.data.frame(v)
f <- function(x) {length(table(cutree(dend, h=x)))}
dm$groupSize <- apply(dm, 1, f)

# Plot dendrogram
dend2 <- hclust(dm_dist, method = clustering_measure)%>% 
  as.dendrogram %>% 
  set_labels(rep("",length(dm_dist))) %>%
  plot()%>%
  abline(h = v, lty = 2, col= colorRampPalette(brewer.pal(8, "Set2"))(length(v)))%>%
  text(x = ncol(dm2) + 20, y = v, labels = dm$groupSize, cex = 1) %>%
  axis(side = 2, at = v, labels = F)

file_name <- paste0("dendrogram2_", lowerLimit, "_", upperLimit, "_", interval, ".png")
dev.copy(png, filename = file_name, res=400, units = "in", width = 10, height = 7)
dev.off()
cat("Plot saved as", file_name, "\n")

################################################
### users should specify height to evaluate based on previous plot ###

height<-0.75

################################################

dend <- hclust(dm_dist, method = "ward.D2") %>% 
  as.dendrogram %>% 
  set_labels(rep("",length(dm_dist))) %>%
  color_branches(h=height) %>% 
  plot()%>%
  abline(h = height, lty = 2, col="black")

file_name <- paste0("dendrogram3_", clustering_measure, "_h", height,".png")
dev.copy(png, filename = file_name, res=400, units = "in", width = 10, height = 7)
dev.off()
cat("Plot saved as", file_name, "\n")

# Obtain motif families
groupsFile <- paste0("MotifFamily_h", height, "_group")

dend <- hclust(dm_dist, method = clustering_measure)
groups<-cutree(dend, h= height)
motifs <- as.vector(read.csv(motifsFile, header = F))
motifs$group<- groups

for (i in c(1:length(table(groups)))){
  groupX <- motifs$V1[motifs$group == i]
  fileName = paste(groupsFile, i, ".tsv", sep = "")
  write.table(groupX, fileName, sep = "\t", quote= F, col.names = F, row.names = F)
}
