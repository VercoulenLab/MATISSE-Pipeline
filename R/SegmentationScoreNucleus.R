library(raster)
library(dplyr)
library(tiff)
library(reshape2)
library(ggplot2)


rootpath <- "/"

## Score for manually annotated nuclei versus automated segmentation of nuclei
# Path for images of IMC and MATISSE detected nuclei outlines
edgePath <- paste0(rootpath, "MapEdge/")

# Path with user annotation nuclei
maskPath <- paste0(rootpath, "/UserMasks//")

SegmentationMethodSuffixIMC <- "_IMC.tiff"
SegmentationMethodSuffixMATISSE <- "_MATISSE.tiff"

regionlist <- list.files(maskPath, full.names = T)

for (i in seq_along(regionlist)) {
  filename <- basename(regionlist[i])
  name <- substr(filename, 0, regexpr(".tiff", filename)-1)

  # Read manually annotated nuclei to binary matrix
  mask <- readTIFF(source = regionlist[i], as.is = T)/255

  # Read automated segmentation masks
  edgeIMC <- readTIFF(paste0(edgePath, name, SegmentationMethodSuffixIMC))
  edgeMATISSE <- readTIFF(paste0(edgePath, name, SegmentationMethodSuffixMATISSE))

  # Calculate probability
  pNuclIMC <- 1 - 1/sum(edgeIMC) * sum(mask * edgeIMC)
  pNuclMATISSE <- 1 - 1/sum(edgeMATISSE) * sum(mask * edgeMATISSE)

  if (i == 1) {
    DF <- list()
  }
  DF[[i]] <- data.frame(i = i, Method = c("IMC", "MATISSE"), Score = c(pNuclIMC, pNuclMATISSE))

  if (i == length(regionlist)) {
    DF <- bind_rows(DF) %>%
      mutate(Method = factor(Method, levels = c("IMC", "MATISSE")))
    write.csv(DF, paste0(rootpath, "DF.csv"))
  }
}

# Create figure for results
png(file = paste0(rootpath, "NuclearCrossScore.png"),
    width = 50, height = 50, units = "mm", res = 1200)
ggplot(DF, aes(x = Method, y = Score))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = .1)+
  geom_line(aes(group = i, color = as.factor(i)), size = .25)+
  theme_bw(base_size = 10)+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        legend.key.size = unit(.5,"line"),
        axis.ticks.length=unit(0, "null"))
dev.off()

# Transform table to wide for t-test
Values <- DF %>%
  tidyr::spread(key = "Method", value = "Score")

# paired T-test
t.test(Values$IMC, Values$MATISSE, paired = T, conf.level = 0.95)
