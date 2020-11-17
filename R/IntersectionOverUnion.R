
### This function is used to calculate Intersection over Union Score (IOU) 
### and Split Merge Analysis.



#Function Argument

IOUCalculation <- function(SegmentationMap_Polygons, Method, imageNumber, scale) {
  SegmentationMap_Polygons_Method <- SegmentationMap_Polygons %>%
    filter(image == imageNumber) %>%
    st_as_sf() %>%
    split(f = .$Method)
  
  IntersectManual_Method <- st_intersects(SegmentationMap_Polygons_Method[["Manual"]], SegmentationMap_Polygons_Method[[Method]], sparse = T)
  
  for (ManualCellIndex in seq_along(IntersectManual_Method)) {
    if (ManualCellIndex == 1) {
      InteractionDF <- list()
      x <- 1
    }
    
    InteractionIndexes <- IntersectManual_Method[[ManualCellIndex]]
    
    ROInrManual <- SegmentationMap_Polygons_Method[["Manual"]][ManualCellIndex, "ROInr"] %>%
      st_drop_geometry() %>%
      as.numeric()
    
    for (InteractionIndex in seq_along(InteractionIndexes)) {
      MethodCellIndex <- InteractionIndexes[InteractionIndex]
      
      ROInrMethod <- SegmentationMap_Polygons_Method[[Method]][MethodCellIndex, "ROInr"] %>%
        st_drop_geometry() %>%
        as.numeric()
      
      ManualPolygon <- SegmentationMap_Polygons_Method[["Manual"]][ManualCellIndex,] %>%
        st_make_valid()
      MethodPolygon <- SegmentationMap_Polygons_Method[[Method]][MethodCellIndex,]
      
      PolyIntersection <- st_intersection(ManualPolygon,
                                          MethodPolygon)
      
      PolyUnion <- st_union(ManualPolygon,
                            MethodPolygon)
      
      PolyIntersectionArea <- st_area(PolyIntersection) / (scale^2)
      PolyUnionArea <- st_area(PolyUnion) / (scale^2)
      IOU <- PolyIntersectionArea / PolyUnionArea
      
      AreaManual <- st_area(ManualPolygon) / (scale^2)
      AreaMethod <- st_area(MethodPolygon) / (scale^2)
      
      InteractionDF[[x]] <- data.frame(image = imageNumber,
                                       ROInrManual,
                                       ROInrMethod,
                                       AreaManual,
                                       AreaMethod,
                                       ManualCellIndex,
                                       MethodCellIndex,
                                       Method = Method,
                                       #PolyIntersection,
                                       PolyIntersectionArea,
                                       #PolyUnion,
                                       PolyUnionArea,
                                       IOU)
      x <- x + 1
    }
    if (ManualCellIndex == length(IntersectManual_Method)) {
      InteractionDF <- bind_rows(InteractionDF)
    }
  }
  return(InteractionDF)
}

# Calculate IOU
Interactions <- IOUCalculation(###)

# Combine IOU data with original polygons of manual annotations
Interactions.MATISSE <- left_join(SegmentationMap_Polygons %>%
                                    filter(Method == "Manual"),
                                  Interactions %>%
                                    filter(Method == "MATISSE"),
                                  by = c("image", "ROInr" = "ROInrManual"))

# Recall per image
Recall_MATISSE <- Interactions.MATISSE %>%
  filter(between(AreaManual, 35, 300)) %>%
  group_by(image, ROInr) %>%
  top_n(n = 1, wt = IOU) %>%
  ungroup() %>%
  group_by(image) %>%
  summarise(TP100 = sum(IOU >= 1.0, na.rm = T),
            TP95 = sum(IOU >= 0.95, na.rm = T),
            TP90 = sum(IOU >= 0.90, na.rm = T),
            TP85 = sum(IOU >= 0.85, na.rm = T),
            TP80 = sum(IOU >= 0.80, na.rm = T),
            TP75 = sum(IOU >= 0.75, na.rm = T),
            TP70 = sum(IOU >= 0.70, na.rm = T),
            TP65 = sum(IOU >= 0.65, na.rm = T),
            TP60 = sum(IOU >= 0.60, na.rm = T),
            TP55 = sum(IOU >= 0.55, na.rm = T),
            TP50 = sum(IOU >= 0.50, na.rm = T),
            TP45 = sum(IOU >= 0.45, na.rm = T),
            TP40 = sum(IOU >= 0.40, na.rm = T),
            TP35 = sum(IOU >= 0.35, na.rm = T),
            TP30 = sum(IOU >= 0.30, na.rm = T),
            Count = n_distinct(ManualCellIndex)) %>%
  gather(key = TP, value = TP_MATISSE, TP100:TP30) %>%
  mutate(Recall _MATISSE = (2 * TP_MATISSE) / (TP_MATISSE + Count),
         TP = as.numeric(substr(TP, 3, nchar(TP)))/100)

# Split per image
MATISSE.Split <- Interactions.MATISSE %>%
  mutate(Overlap = PolyIntersectionArea / AreaManual) %>%
  filter(Overlap >= .2) %>%
  group_by(image, ROInr) %>%
  count() %>%
  ungroup() %>%
  group_by(image) %>%
  summarise(mean = mean(n), median = median(n), max = max(n), objectcount = n(), fraction = sum(n > 1) / objectcount)
