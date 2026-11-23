library(Seurat)
library(spacexr)
library(tidyverse)
library(data.table)
library(mltools)

normalize <- function(x){ ## ordinary 0 to 1 normalization
  x <- (x - base::min(x)) / (max(x) - base::min(x))
  return(x)
}

## D4 ##
rctd.res <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/RCTD/day4/Day4mouseRCTDoutmin30_editv1")
## D10 ##
rctd.res <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/RCTD/day10/Day10mouseRCTDoutmin30_editv1")
## D4 full mode##
rctd.res <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/Day4mouseRCTDoutmin30_multi.rds")
## D10 full mode##
rctd.res <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD_day10/Day10mouseRCTDoutmin30_multi.rds")

rctd.res <- rctd.res@RCTD.reps

res <- list()

for(i in 1:2){
  
  ## get cell type weights
  map = as.data.frame(as.matrix(rctd.res[[i]]@results$weights))
  
  ## format
  map = as.data.frame(apply(X = map, MARGIN = 2, FUN = normalize))
  
  ## get coords
  coords = as.data.frame(rctd.res[[i]]@spatialRNA@coords)
  coords = coords[rownames(map),]
  
  ## centre coords for rotation
  coords$x = coords$x - mean(coords$x)
  coords$y = coords$y - mean(coords$y)
  

  
  ## get coord data
  map$x = coords$x
  map$y = coords$y
  map = map %>% relocate(y) %>% relocate(x)
  
  ## write
  res[[i]] <- map
}

## plot results
ggplot() +
  geom_point(
    data = res[[1]], mapping = aes(x=x,y=y,alpha=`T cells`),
    color = "green", size = 0.5) +
  geom_point(
    data = res[[1]], mapping = aes(x=x,y=y,alpha=`B cells`),
    color = "red", size = 0.5) +
  geom_point(
    data = res[[1]], mapping = aes(x=x,y=y,alpha=`Epithelial Cells`),
    color = "cyan", size = 0.5) +
  theme_minimal() +
  scale_alpha_continuous(range = c(0,1), breaks=c(0,0.1,1)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    axis.line = element_blank(),
    text = element_blank(),
    plot.title = element_text(size = 30),
    plot.title.position = 'panel',
    axis.title = element_blank()) +
  #xlim(c(600, 1300)) + ylim(c(-400,300)) +
  coord_fixed()


