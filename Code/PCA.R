own_ss <- readRDS("rds/own_sample_sheet.rds")
data <- readRDS("rds/Freichel.rds")
library(DESeq2)

#x <- own_ss$data_Names
#own_ss[1] <- rownames(own_ss)
#rownames(own_ss) <- x

#ddsMat <- DESeqDataSetFromMatrix(countData = data[2:97], colData = own_ss,design = ~1 )

#rld <- rlog(ddsMat)

##### part above only needed if rld.rds has to be computed

rld <- readRDS("rds/rdl.rds")

png('./Figures/PCA/PCA.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld,"Info_treatment_group")

dev.off()


############################## adding columns for seperation of groups...
colData(rld)$Genotype <- "NA"
colData(rld)$Genotype[grepl("Stim1",colData(rld)$Info_treatment_group)] <- "Stim1/2_control"
colData(rld)$Genotype[grepl("Wu4",colData(rld)$Info_treatment_group)] <- "Wu4_control"
colData(rld)$Genotype[grepl("WT",colData(rld)$Info_treatment_group)] <- "WT"
colData(rld)$Genotype[grepl("Stim1.+KO",colData(rld)$Info_treatment_group)] <- "Stim1/2_KO"
colData(rld)$Genotype[grepl("Wu4.+KO",colData(rld)$Info_treatment_group)] <- "Wu4_KO"
colData(rld)$Genotype[grepl("TRPC1.+KO",colData(rld)$Info_treatment_group)] <- "Wu4_KO"



colData(rld)$Treatment <- "NA"
colData(rld)$Treatment[grepl("saline",colData(rld)$Info_treatment_group)] <- "saline"
colData(rld)$Treatment[grepl("Isoproterenol",colData(rld)$Info_treatment_group)] <- "Isoproterenol"
colData(rld)$Treatment[grepl("Angiotensin",colData(rld)$Info_treatment_group)] <- "Angiotensin"
colData(rld)$Treatment[grepl("sham",colData(rld)$Info_treatment_group)] <- "sham"
colData(rld)$Treatment[grepl("TAC",colData(rld)$Info_treatment_group)] <- "TAC_48_hours"
colData(rld)$Treatment[grepl("TAC.+day",colData(rld)$Info_treatment_group)] <- "TAC_(7_days)"



colData(rld)$p_t <- "NA"
colData(rld)$p_t[grepl("p",colData(rld)$V3)] <- "p"
colData(rld)$p_t[grepl(" t",colData(rld)$V3)] <- "t"







##############################




### Subset plots


###Genotype
png('./Figures/PCA/PCA_Genotype.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld,"Genotype")

dev.off()

###Treatment
png('./Figures/PCA/PCA_Treatment.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld,"Treatment")

dev.off()

###p/t
png('./Figures/PCA/PCA_p_t.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld,"p_t")

dev.off()



#Stim1 Stim2
o <- grepl("Stim",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Stim.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Stim1 Stim2 saline
o <- grepl("Stim.+saline",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Stim_saline.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Stim1 Stim2 Angiotensin
o <- grepl("Stim.+Ang",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Stim_ATII.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()



#TRPC KO / WT
o <- grepl("TRPC|WT",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_TRPC.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#TRPC sham
o <- grepl("TRPC.+sham|WT.+sham",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_TRPC_sham.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#TRPC sham
o <- grepl("TRPC.+TAC|WT.+TAC",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_TRPC_TAC.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Wu4
o <- grepl("Wu4",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Wu4.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()


#Wu4 Saline
o <- grepl("Wu4.+saline",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Wu4_saline.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Wu4 Isoproterenol
o <- grepl("Wu4.+Iso",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Wu4_Iso.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Wu4 Sham
o <- grepl("Wu4.+sham",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Wu4_sham.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()

#Wu4 TAC
o <- grepl("Wu4.+TAC",own_ss$Info_treatment_group)

png('./Figures/PCA/PCA_Wu4_TAC.png',width = 40, height = 20, units = "cm", res=90)

plotPCA(rld[,o],"Info_treatment_group")

dev.off()


#additional groups

#png('./Figures/PCA/PCA_Control_KO.png',width = 40, height = 20, units = "cm", res=90)
#colData(rld)$Genotype <- "Treat"
#colData(rld)$Genotype[grepl("Genotype",colData(rld)$Info_treatment_group)] <- "Genotype"

#plotPCA(rld,"Genotype")

#dev.off()


