# creating DESeq2 dataset with appropriate design and MA plot the data

own_ss <- readRDS("rds/own_sample_sheet.rds")
data <- readRDS("rds/Freichel.rds")
library(DESeq2)
library(openxlsx)

# loading and preparing gene alignement
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
fck <- getBM(attributes = c("ensembl_gene_id","description"), filters="ensembl_gene_id",values=data$gene_id,mart = ensembl)
colnames(fck)[1] <- "gene_id"

#renaming columns
x <- own_ss$data_Names
own_ss[1] <- rownames(own_ss)
rownames(own_ss) <- x

#only the Wu4 subset is used

o <- grepl("Wu4",own_ss$Genotype)
u <- grepl("p",own_ss$p_t)
m <- u & o
helper <-  logical(1) 
helper[2:97] <- m

own_ss$Treatment <- relevel(factor(own_ss$Treatment),"saline")

deseqMat <-
  DESeqDataSetFromMatrix(
    countData = data[helper],
    colData = own_ss[m, ],
    design = ~ Genotype + Treatment + Genotype:Treatment
  )
ppp <- DESeq(deseqMat)

res <- results(ppp)

plotMA(res, main="DESeq2", ylim=c(-2,2))

##############
#the usage of the interaction terms needs to be checked! 
#resultsNames(ppp)

#################


# Control Treatment_saline_isoprot effect
res <- results(ppp, contrast = c("Treatment", "saline", "Isoproterenol"))

png('./Figures/MA_Plot/Wu4/Control_Treatment_saline_isoprot.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/Control_Treatment_saline_isoprot.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



# Control Treatment_saline_sham effect
res <- results(ppp, contrast = c("Treatment", "saline", "sham"))

png('./Figures/MA_Plot/Wu4/Control_Treatment_saline_sham.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/Control_Treatment_saline_sham.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



# Control Treatment_saline_Tac(48) effect
res <- results(ppp, contrast = c("Treatment", "saline", "TAC_48_hours"))

png('./Figures/MA_Plot/Wu4/Control_Treatment_saline_Tac(48).png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/Treatment_saline_Tac(48).xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")




# Control Treatment_saline_Tac(7d) effect
res <- results(ppp, contrast = c("Treatment", "saline", "TAC_(7_days)"))

png('./Figures/MA_Plot/Wu4/Control_Treatment_saline_Tac(7d).png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/Treatment_saline_Tac(7d).xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



#KO Treatment_saline_isoprot effect
res <- results(ppp,name = "GenotypeWu4_KO.TreatmentIsoproterenol")

png('./Figures/MA_Plot/Wu4/KO_Treatment_saline_isoprot.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/KO_Treatment_saline_isoprot.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



#KO Treatment_saline_sham effect ??????????????????? saline vs sham??
res <- results(ppp,name = "GenotypeWu4_KO.Treatmentsham")

png('./Figures/MA_Plot/Wu4/KO_Treatment_saline_sham.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/KO_Treatment_saline_sham.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



#KO Treatment_saline_TAC(48h) effect
res <- results(ppp,name = "GenotypeWu4_KO.TreatmentTAC_48_hours")

png('./Figures/MA_Plot/Wu4/KO_Treatment_saline_TAC(48h).png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/KO_Treatment_saline_TAC(48h).xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



#KO Treatment_saline_TAC(7d) effect
res <- results(ppp,name = "GenotypeWu4_KO.TreatmentTAC_.7_days.")

png('./Figures/MA_Plot/Wu4/KO_Treatment_saline_TAC(7d).png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/KO_Treatment_saline_TAC(7d).xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")




#difference between treatment effect on genotypes
res <- results(ppp, name = "Genotype_Wu4_KO_vs_Wu4_control")
png('./Figures/MA_Plot/Wu4/Genotype_treatmenteffectdiff.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")#, ylim=c(-2,2))
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Wu4/Genotype_Treatmenteffectdiff.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")


  
