#PCA-Plot with ggplot

library(ggplot2)
own_ss <- readRDS("rds/own_sample_sheet.rds")
data <- readRDS("rds/Freichel.rds")

q <- prcomp(t(log2(1 + (data[2:97]))))
df <- as.data.frame(q$x)

ggplot(df,aes(x=df$PC1,y=df$PC2, colour = own_ss$Genotype))  + geom_point() + xlab("PC1") + ylab("PC2") 

#+geom_text(aes(label=Name),hjust=0, vjust=0)


###########################################################
# Plot only saline with names 
o <- grepl("saline",own_ss$Treatment)
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = own_ss[o,]$Treatment
))  + geom_point() + xlab("PC1") + ylab("PC2")  +geom_text(aes(label= own_ss$V3[o]),hjust=0, vjust=0)

#Plot without t
o <- !(grepl("t",own_ss$p_t))
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = own_ss[o,]$p_t
))  + geom_point() + xlab("PC1") + ylab("PC2")  +geom_text(aes(label= own_ss$V3[o]),hjust=0, vjust=0)

#Plot Genotype without t
o <- !(grepl("t",own_ss$p_t))
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = own_ss[o,]$Genotype
))  + geom_point() + xlab("PC1") + ylab("PC2")








#######################################################
#Plot Genotype without t
Genotype <- own_ss[o,]$Genotype

png('./Figures/PCA/PCA_Genotype_wo_t.png',width = 40, height = 20, units = "cm", res=90)

o <- !(grepl("t",own_ss$p_t))
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = Genotype
))  + geom_point() + xlab("PC1") + ylab("PC2")

dev.off()

#Plot Treatment without t
Treatment <- own_ss[o,]$Treatment
png('./Figures/PCA/PCA_Treatment_wo_t.png',width = 40, height = 20, units = "cm", res=90)

o <- !(grepl("t",own_ss$p_t))
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = Treatment
))  + geom_point() + xlab("PC1") + ylab("PC2")


dev.off()

#Show the outlier
Treatment <- own_ss[o,]$Treatment
png('./Figures/PCA/PCA_Treatment_wo_t.png',width = 40, height = 20, units = "cm", res=90)

w <- grepl("saline",own_ss$Treatment)
h <- (grepl("t",own_ss$p_t))
o <- h & (grepl("Wu4_control",own_ss$Genotype)) & w
ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = own_ss[o,]$Genotype
))  + geom_point() + xlab("PC1") + ylab("PC2") +geom_text(aes(label=own_ss[o,]$V3),hjust=0, vjust=0)
