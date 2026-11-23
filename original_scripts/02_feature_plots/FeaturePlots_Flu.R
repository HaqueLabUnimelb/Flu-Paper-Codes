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

## Megre objects, for VlnPlot purposes

mouse1$sample <- "mouse 4.1"
mouse2$sample <- "mouse 4.2"
mouse4$sample <- "mouse 10.4"
mouse5$sample <- "mouse 10.5"

combined <- merge(mouse1,
                  y = list(mouse2, mouse4, mouse5),
                  add.cell.ids = c("mouse 4.1","mouse 4.2","mouse 10.4","mouse 10.5"))

combined$sample <- factor(combined$sample,
                          levels = c("mouse 4.1","mouse 4.2","mouse 10.4","mouse 10.5"))

DefaultAssay(combined) <- "RNA"

library(dplyr)

expr <- FetchData(combined, vars = c("Il1b", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Il1b > 0) * 100)

expr <- FetchData(combined, vars = c("Il6", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Il6 > 0) * 100)

expr <- FetchData(combined, vars = c("Cxcl1", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Cxcl1 > 0) * 100)

expr <- FetchData(combined, vars = c("Tnf", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Tnf > 0) * 100)

expr <- FetchData(combined, vars = c("Il10", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Il10 > 0) * 100)

expr <- FetchData(combined, vars = c("Ifnb1", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Ifnb1 > 0) * 100)

expr <- FetchData(combined, vars = c("Ifng", "sample"))

expr %>%
  group_by(sample) %>%
  summarise(percent_detected = mean(Ifng > 0) * 100)


VlnPlot(combined, features = "Il1b", group.by = "sample", pt.size = 1)+NoLegend() 
ggsave(filename = "Il1b_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Il6", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Il6_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Cxcl1", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Cxcl1_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Tnf", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Tnf_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Ifng", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Ifng_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Ifnb1", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Ifnb1_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")

VlnPlot(combined, features = "Il10", group.by = "sample", pt.size = 1)+NoLegend()
ggsave(filename = "Il10_VlnPlot.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/vlnplots/", width = 12, height = 12, units = "cm")


VlnPlot(combined, features = "nFeature_RNA", group.by = "sample", pt.size = 0)


## Plot
FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'nFeature_RNA', order = TRUE, raster = FALSE) +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(
    option = 'magma',
    direction = -1,
    begin = 0.1,
    end = 0.95
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

FeaturePlot(mouse1, reduction = 'SPATIAL', features = 'Il1b', order=T, raster = F) + coord_fixed() + 
  theme(axis.ticks=element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = "white")) +
  scale_colour_viridis_c(option = 'magma', direction = -1, begin = 0.1, end = 0.95)+ NoAxes()
ggsave(filename = "mouse4_ifnb1_v2.pdf", device = "pdf", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplots/", width = 12, height = 12, units = "cm")

mouse1$gene_plot <- FetchData(mouse1, vars = "Il1b")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Il6")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Cxcl1")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Tnf")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Il10")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Ifnb1")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse1$gene_plot <- FetchData(mouse1, vars = "Ifng")[,1]
mouse1$gene_plot[mouse1$gene_plot == 0] <- NA

FeaturePlot(mouse1, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse1_ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")



### mouse2

mouse2$gene_plot <- FetchData(mouse2, vars = "Il1b")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Il6")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Cxcl1")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Tnf")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Il10")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Ifnb1")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse2$gene_plot <- FetchData(mouse2, vars = "Ifng")[,1]
mouse2$gene_plot[mouse2$gene_plot == 0] <- NA

FeaturePlot(mouse2, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse2_ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")


### mouse4

mouse4$gene_plot <- FetchData(mouse4, vars = "Il1b")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Il6")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Cxcl1")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Tnf")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Il10")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Ifnb1")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse4$gene_plot <- FetchData(mouse4, vars = "Ifng")[,1]
mouse4$gene_plot[mouse4$gene_plot == 0] <- NA

FeaturePlot(mouse4, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse4_ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")


### mouse5

mouse5$gene_plot <- FetchData(mouse5, vars = "Il1b")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_il1b.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Il6")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_il6.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Cxcl1")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_cxcl1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Tnf")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_tnf.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Il10")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_il10.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Ifnb1")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_ifnb1.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")

##
mouse5$gene_plot <- FetchData(mouse5, vars = "Ifng")[,1]
mouse5$gene_plot[mouse5$gene_plot == 0] <- NA

FeaturePlot(mouse5, reduction = "SPATIAL", features = "gene_plot",
            order = TRUE, raster = FALSE) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1,
    begin = 0.1,
    end = 0.95,
    na.value = "grey90"
  ) +
  NoAxes()+ coord_fixed()

ggsave(filename = "mouse5_ifng.jpeg", device = "jpeg", path = "/data/projects/punim2394/lung_curio_mouse_pr8/revision_analysis/featureplotsv2/", width = 10, height = 10, units = "cm")
