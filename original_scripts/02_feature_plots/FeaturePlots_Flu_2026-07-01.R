library(Seurat)
library(spacexr)
library(tidyverse)
library(data.table)
library(mltools)

# Day 4, mouse 1
mouse1 <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_2_2_QCed.rds")
# Day 4, mouse 2
mouse2 <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_2_3_QCed.rds")
# Day 10, mouse 4
mouse4 <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_4_2_QCed.rds")
# Day 10, mouse 5
mouse5 <- readRDS("/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_4_3_QCed.rds")

mouse1 <- NormalizeData(mouse1)
mouse2 <- NormalizeData(mouse2)
mouse4 <- NormalizeData(mouse4)
mouse5 <- NormalizeData(mouse5)

## plot



FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'nFeature_RNA', order = TRUE, raster = FALSE) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(
    option = 'magma',
    direction = -1,
    begin = 0.1,
    end = 0.95,
    trans = "log10" 
  )+NoAxes()
ggsave(filename = "mouse1_umis.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/", width = 12, height = 12, units = "cm")


FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'nFeature_RNA', order = TRUE, raster = FALSE) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(
    option = 'magma',
    direction = -1,
    begin = 0.1,
    end = 0.95,
    trans = "log10" 
  )+NoAxes()
ggsave(filename = "mouse2_umis.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/", width = 12, height = 12, units = "cm")


FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'nFeature_RNA', order = TRUE, raster = FALSE) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(
    option = 'magma',
    direction = -1,
    begin = 0.1,
    end = 0.95,
    trans = "log10" 
  )+NoAxes()
ggsave(filename = "mouse4_umis.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'nFeature_RNA', order = TRUE, raster = FALSE) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(
    option = 'magma',
    direction = -1,
    begin = 0.1,
    end = 0.95,
    trans = "log10" 
  )+NoAxes()
ggsave(filename = "mouse5_umis.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/", width = 12, height = 12, units = "cm")


# Gene expression of BAL cytokines in Curioseeker arrays - Tnf, Ifng, Il6, Ifna1, Ifnb1, Cxcl1, Il10, Il1b

# Mouse1
FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Tnf', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Ifng', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Il6', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# could not detect expression of Ifna1 in array
FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Ifna1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Ifnb1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Cxcl1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Il10', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Il1b', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

#Mouse2
FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Tnf', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Ifng', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Il6', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# could not detect expression of Ifna1 in array
FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Ifna1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Ifnb1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Cxcl1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Il10', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Il1b', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

#Mouse4
FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Tnf', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Ifng', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Il6', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# could not detect expression of Ifna1 in array
FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Ifna1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

# could not detect expression of Ifnb1 in array
FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Ifnb1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Cxcl1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Il10', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Il1b', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

#Mouse5
FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Tnf', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Ifng', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Il6', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# could not detect expression of Ifna1 in array
FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Ifna1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

# could not detect expression of Ifnb1 in array
FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Ifnb1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Cxcl1', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Il10', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Il1b', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

#~~~~~~~~~~~~~
#~~~~~~~~~~~~~
#~~~~~~~~~~~~~

# Plot B clelr elated genes - most spatially variable features spat out from Morans I analysis Curioseeker pipeline

# Mouse 1
FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Igkc', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Igkc.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Ighm', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Ighm.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Ighg2c', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse1_Ighg2c.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# Mouse 2
FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Igkc', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Igkc.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Ighm', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse2_Ighm.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# Not featured
FeaturePlot(mouse2, reduction = 'SPATIAL', features = 'Ighg2c', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()

# Mouse 4
FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Igkc', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Igkc.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Ighm', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Ighm.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse4, reduction = 'SPATIAL', features = 'Ighg2c', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_Ighg2c.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

# Mouse 5
FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Igkc', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Igkc.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Ighm', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Ighm.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

FeaturePlot(mouse5, reduction = 'SPATIAL', features = 'Ighg2c', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse5_Ighg2c.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")












