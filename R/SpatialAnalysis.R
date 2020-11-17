######### This Script performs Spatial Analysis. #######


#Packages required to be installed & loaded before analysis 
library(raster)
library(stars)
library(sf)
library(magrittr)
library(dplyr)
library(doParallel)
library(RANN)

projectpath <- "/"

### Import single-cell data
load(paste0(projectpath, "SingleCellData-MATISSE.Rda"))

# For parallel processing define number of aviable cores
#aviablecores <- 8
#registerDoParallel(cores = aviablecores)

getPoly <- function(path){
  polys <- raster(path) %>%
    na_if(0) %>%
    st_as_stars() %>%
    st_as_sf() %>%
    group_by_at(1) %>%
    summarise() %>%
    rename_at(1, ~ "ROInr")
  polys$path <- path
  polys$ImageNumber <- substr(basename(path), 0, regexpr(SegmentationMapSuffix, basename(path))-1)
  return(polys)
}

# Folder where segmentation maps are located
SegmentationPath <- paste0(projectpath, "SegmentationMaps/Objects/")
SegmentationMapSuffix <- "_MATISSE_Cells.tiff"

# List all segmentation maps
Filelist <- data.frame(path = list.files(path = SegmentationPath, pattern = paste0(SegmentationMapSuffix, "$"), full.names = T))
Filelist$ImageNumber <- substr(basename(Filelist$path), 0, regexpr(SegmentationMapSuffix, basename(Filelist$path))-1)

# Keep segmentation map paths only for regions present in single-cell data
Filelist$path[is.element(Filelist$ImageNumber, FusedDF$ImageNumber)]

# Parallel retrieval of polygons for all cells in all segmentation maps
MapPolygons <- foreach(i = seq_along(Filelist.unique), .combine = rbind, .packages = c("dplyr","raster","sf","stars")) %dopar% {
  getPoly(Filelist.unique[i])
}

# Join polygons with single-cell data and save to new file
FusedDF.polys <- inner_join(FusedDF, MapPolygons, by = c("ImageNumber", "ROInr"))
save(FusedDF.polys, file = paste0(projectpath, "SingleCellData-MATISSE_polys.Rda"))


######## Calculate Distance between neighbouring cells #######

# Select crucial columns
TempDF <- FusedDF.polys %>%
  dplyr::select(ImageNumber, ROInr, geometry)

ImageNRS <- as.numeric(unique(TempDF$ImageNumber))

NumberInRange <- function(DF, radius) {
  distances <- DF %>%
    as.data.frame() %>%
    dplyr::select(contains("nn.dists"))

  for (i in seq_along(radius)) {
    DF[, paste0("nr.in.range.", radius[i])] <- rowSums(distances < radius[i])
  }
  return(DF)
}

AddSpatialInfo <- function(DF, radius) {
  DF <- st_as_sf(DF)
  DF[,"centroid"] <- DF %>%
    dplyr::select(geometry) %>%
    st_centroid()
  centroids <- do.call(rbind, DF$centroid)
  closest <- nn2(centroids, k = 26, searchtype = "standard") %>%
    as.data.frame() %>%
    # exclude self-self distance
    dplyr::select(-c(nn.idx.1,nn.dists.1))

  DF <- cbind(DF, closest)
  DF <- NumberInRange(DF, radius)
  return(DF)
}

for (i in seq_along(ImageNRS)) {
  print(paste("Processing image:", i))
  temp <- TempDF %>%
    filter(ImageNumber == ImageNRS[i])
    if (i == 1) {
      TempDF.suppl <- AddSpatialInfo(temp, c(5,10,15,20,25))
    } else {
      TempDF.suppl <- bind_rows(TempDF.suppl, AddSpatialInfo(temp, c(5,10,15,20,25)))
    }
}

FusedDF.polys.distances <- inner_join(FusedDF, TempDF.suppl %>%
                                       as.data.frame() %>%
                                       dplyr::select(- geometry), by = c("ImageNumber", "ROInr"))

save(FusedDF.polys.distances, file = paste0(projectpath, "SingleCellData-MATISSE_polys_distances.Rda"))


