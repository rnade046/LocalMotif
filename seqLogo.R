#Package dependencies
packages = c("ggseqlogo", "ggplot2")

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
groups <- args[3]
height <- args[4]

# wd <- "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\corrNetTop2-400_coreTPD_p0.4_p5.41545109270352E-7\\Groups_h0.9\\"
# projectName <- "corrNetTop2-400_coreTPD_p0.4"
# groups <- 18
# height <- 0.9
setwd(wd)

data <- list()
for(i in 1:groups){
  
  file <- paste(projectName, "_h", height,"_motifInstances_motifFamily", i, sep = "")
  motifs <- read.csv(file, sep = "\t", header = F)
  
  data[[length(data)+1]] <- unlist(motifs)
  
}

outputFile <- paste(projectName, "_h", height, "_MotifFamilies.png")

names(data) <- seq(1, groups, 1)

# ggseqlogo(data, ncol = 2, method = "prob") + 
#   theme(axis.text.x=element_text(colour="white"), 
#         axis.text.y = element_text(color = "white"),
#         text = element_text(colour = "white"),
#         strip.text.x = element_text(colour = "white"))
# ggsave(outputFile, width = 5, height = 10,  bg = "transparent")

ggseqlogo(data, ncol = 2, method = "prob")
ggsave(outputFile, width = 5, height = 10,  bg = "transparent")
