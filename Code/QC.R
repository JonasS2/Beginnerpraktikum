own_ss <- readRDS("rds/own_sample_sheet.rds")
data <- readRDS("rds/Freichel.rds")

colnames(data)[2:97] <-own_ss$V3 # change to useful colnames in original dataset


sampleDists <- dist(log2(1 + t(data[2:97])))

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


png('./Figures/heatmap.png',width = 30, height = 30, units = "cm", res=90)

pheatmap(sampleDistMatrix,
              clustering_distance_rows = sampleDists,
              clustering_distance_cols = sampleDists,
              col = colors)

dev.off()